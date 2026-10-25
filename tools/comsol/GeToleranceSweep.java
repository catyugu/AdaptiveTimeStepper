import com.comsol.model.*;
import com.comsol.model.util.*;
import java.io.*;
import java.util.*;

// Which axes of COMSOL's absolute tolerance move the trajectory, and how?
//
// Two models, one per group of cases:
//
//   field  A 1D Coefficient Form PDE on the default unit interval with five
//          linear elements plus one Global Equation, so the model has two
//          variables of very different size:
//            u: da=1, a=1, c=0, f=0, u(x,0) = 1+x  ->  dofs 1.0 ... 2.0
//            v: Global Equations, vt+0.5*v = 0, v(0) = 0.01  ->  1 dof
//          c=0 and a=1 give M*u' + M*u = 0 with the same mass matrix on both
//          sides, so every dof obeys u_i(t) = u_i(0)*exp(-(t-t0)) exactly and
//          the exported state has an analytic reference:
//            u(x,t) = (1+x)*exp(-(t-t0)),  v(t) = 0.01*exp(-0.5*(t-t0)).
//          The initial value is not uniform on purpose: the per-dof magnitudes
//          differ, so a per-dof weight, a per-variable weight and a manual
//          absolute tolerance are three different numbers.
//
//   decay  The single Global Equation of comsol_twofield.hpp (u'+u=0, u(0)=1,
//          span 10). These cases exist to reconcile the two recorded first
//          steps, which disagree because one set pins the variable scale and
//          the other leaves it automatic.
//
// Every case writes name.csv (time, all dofs, all dof derivatives), name.log
// (the solver log, which carries the scale block and the accepted-step table)
// and one manifest.csv row holding every value read back from the solver.
//
// The manual-scale, non-unit cases are the ones that separate the two readings
// of the scaled tolerance: whether the variable's scale divides the solution in
// the weight (Y = U/sigma) or only multiplies the absolute tolerance. With
// sigma = 1 everywhere those readings coincide, which is why every earlier
// sweep agreed with both.
public class GeToleranceSweep {
    static final String[] SPECS = {
        // The six cases the sweep was asked for. `method` is atolglobalmethod,
        // `value` is atolglobalvaluemethod, `rescale` is ewtrescale.
        "field_scaled_factor_on",
        "field_scaled_factor_noupdate rescale=off",
        "field_scaled_manual1e4 value=manual atol=1e-4",
        "field_unscaled_factor method=unscaled",
        "field_unscaled_manual1e4 method=unscaled value=manual atol=1e-4",
        // The span question: the same tolerance settings with the start at t=2
        // and rtol=1, so the derivative guard does not clip the initial step and
        // the 0.1% cap itself is the observable.
        "field_span_probe rtol=1 t0=2 end=12",
        // Scale probes. Pinning the variable scale to 1 takes the automatic
        // scale out of the picture: that is the setting the two-field evidence
        // used, and it is the setting the recorded trajectories did not use.
        "field_scale1_on scale=manual",
        "field_scale1_off scale=manual rescale=off",
        "field_scale1_unscaled scale=manual method=unscaled",
        // The scale itself, at values an exactly representable double holds and
        // the readback returns in full. The scaled criterion divides the
        // solution by it; the unscaled criterion is documented not to, so the
        // unscaled pair must reproduce the scale-1 case exactly.
        "field_manualscale2 scale=manual scaleu=2 scalev=0.5",
        "field_manualscale2_unscaled scale=manual scaleu=2 scalev=0.5 method=unscaled",
        // Determinism control: a second solve of the first case.
        "field_scaled_factor_on_repeat",
        // The decay model. `decay_default` is comsol_stepping.hpp "auto_span10";
        // `decay_scale1` is comsol_twofield.hpp "single_a" with its manual scale.
        "decay_default model=decay order=1",
        "decay_scale1 model=decay order=1 scale=manual",
        "decay_scale1_off model=decay order=1 scale=manual rescale=off",
        "decay_auto_off model=decay order=1 rescale=off",
        "decay_scale1_manualatol model=decay order=1 scale=manual value=manual atol=1e-3",
    };

    static String directory;
    static PrintWriter manifest;

    static String name;
    static String model;
    static double t0;
    static double end;
    static String method;
    static String valueMethod;
    static double atolValue;
    static String rescale;
    static String scaleMethod;
    static double scaleU = 1.0;
    static double scaleV = 1.0;
    static double rtol;
    static int maxOrder;
    static int elements;

    static void reset() {
        name = "case";
        model = "field";
        t0 = 0.0;
        end = 10.0;
        method = "scaled";
        valueMethod = "factor";
        atolValue = 1e-3;
        rescale = "on";
        scaleMethod = "auto";
        scaleU = 1.0;
        scaleV = 1.0;
        rtol = 0.01;
        maxOrder = 5;
        elements = 5;
    }

    // "name key=value ..."; every omitted key keeps the default above, so a row
    // shows only what differs from a plain transient solve.
    static void apply(String spec) {
        reset();
        String[] tokens = spec.trim().split("\\s+");
        name = tokens[0];
        for (int i = 1; i < tokens.length; ++i) {
            int split = tokens[i].indexOf('=');
            if (split < 0) throw new IllegalArgumentException("bad token: " + tokens[i]);
            String key = tokens[i].substring(0, split);
            String value = tokens[i].substring(split + 1);
            if (key.equals("model")) {
                model = value;
            } else if (key.equals("t0")) {
                t0 = Double.parseDouble(value);
            } else if (key.equals("end")) {
                end = Double.parseDouble(value);
            } else if (key.equals("method")) {
                method = value;
            } else if (key.equals("value")) {
                valueMethod = value;
            } else if (key.equals("atol")) {
                atolValue = Double.parseDouble(value);
            } else if (key.equals("rescale")) {
                rescale = value;
            } else if (key.equals("scale")) {
                scaleMethod = value;
            } else if (key.equals("scaleu")) {
                scaleU = Double.parseDouble(value);
            } else if (key.equals("scalev")) {
                scaleV = Double.parseDouble(value);
            } else if (key.equals("rtol")) {
                rtol = Double.parseDouble(value);
            } else if (key.equals("order")) {
                maxOrder = Integer.parseInt(value);
            } else if (key.equals("elements")) {
                elements = Integer.parseInt(value);
            } else {
                throw new IllegalArgumentException("unknown key: " + key);
            }
        }
        if (!model.equals("field") && !model.equals("decay")) {
            throw new IllegalArgumentException("unknown model: " + name);
        }
        if (!method.equals("scaled") && !method.equals("unscaled")) {
            throw new IllegalArgumentException("unknown atol method: " + name);
        }
        if (!valueMethod.equals("factor") && !valueMethod.equals("manual")) {
            throw new IllegalArgumentException("unknown atol value method: " + name);
        }
        if (!rescale.equals("on") && !rescale.equals("off")) {
            throw new IllegalArgumentException("unknown rescale setting: " + name);
        }
        if (!scaleMethod.equals("auto") && !scaleMethod.equals("manual")) {
            throw new IllegalArgumentException("unknown scale method: " + name);
        }
        if (scaleMethod.equals("auto") && (scaleU != 1.0 || scaleV != 1.0)) {
            throw new IllegalArgumentException("a specified scale needs scale=manual: " + name);
        }
        if (maxOrder < 1 || maxOrder > 5) {
            throw new IllegalArgumentException("order must be 1-5: " + name);
        }
    }

    static String num(double value) {
        if (value == Math.rint(value)) return Long.toString((long) value);
        return Double.toString(value);
    }

    static String csv(String value) {
        return "\"" + value.replace("\"", "\"\"") + "\"";
    }

    static Model build() {
        double step = (end - t0) / 10.0;
        String times = "range(" + num(t0) + "," + num(step) + "," + num(end) + ")";
        Model m = ModelUtil.create("Model");
        m.component().create("comp1", true);
        m.component("comp1").geom().create("geom1", 1);
        m.component("comp1").geom("geom1").create("i1", "Interval");
        m.component("comp1").geom("geom1").run();
        m.component("comp1").mesh().create("mesh1");
        if (model.equals("field")) {
            // hmax = 1/elements gives exactly `elements` elements on [0,1], which
            // the run verifies from getNumElem().
            m.component("comp1").mesh("mesh1").feature("size").set("custom", "on");
            m.component("comp1")
                .mesh("mesh1").feature("size").set("hmax", Double.toString(1.0 / elements));
        }
        m.component("comp1").mesh("mesh1").run();
        m.component("comp1").physics().create("ge", "GlobalEquations", "geom1");
        if (model.equals("field")) {
            m.component("comp1").physics("ge").feature("ge1").set("equation", "vt+0.5*v");
            m.component("comp1").physics("ge").feature("ge1").set("name", "v");
            m.component("comp1").physics("ge").feature("ge1").set("initialValueU", "0.01");
            m.component("comp1").physics().create("c", "CoefficientFormPDE", "geom1");
            m.component("comp1").physics("c").feature("cfeq1").set("da", "1");
            m.component("comp1").physics("c").feature("cfeq1").set("a", "1");
            m.component("comp1").physics("c").feature("cfeq1").set("c", "0");
            m.component("comp1").physics("c").feature("cfeq1").set("f", "0");
            // The non-uniform initial value: the dofs are 1, 1.2, ... 2.
            m.component("comp1").physics("c").feature("init1").set("u", "1+x");
            m.component("comp1").physics("c").prop("ShapeProperty").set("order", "1");
        } else {
            m.component("comp1").physics("ge").feature("ge1").set("equation", "ut+u");
            m.component("comp1").physics("ge").feature("ge1").set("name", "u");
            m.component("comp1").physics("ge").feature("ge1").set("initialValueU", "1");
        }
        m.study().create("std1");
        m.study("std1").create("time", "Transient");
        m.study("std1").feature("time").set("tlist", times);
        m.study("std1").createAutoSequences("all");
        return m;
    }

    // The solver's variable scale nodes, by the role each variable plays: the
    // Coefficient Form PDE's field is comp1_u and a Global Equations variable is
    // comp1_ODE1, whatever it is called.
    static String[] variables() {
        if (model.equals("field")) return new String[]{"comp1_u", "comp1_ODE1"};
        return new String[]{"comp1_ODE1"};
    }

    // The solver's own readback of the variable scales: "manual:1" pins them,
    // "parent:..." means the solver computes them (and prints them in the log).
    static String readScale(Model m, String variable) {
        SolverFeature feature = m.sol("sol1").feature("v1").feature(variable);
        return feature.getString("scalemethod") + ":" + feature.getString("scaleval");
    }

    static void run() throws Exception {
        Model m = build();
        int meshed = model.equals("field") ? m.component("comp1").mesh("mesh1").getNumElem() : 0;
        if (model.equals("field") && meshed != elements) {
            throw new IllegalStateException("expected " + elements + " elements, meshed "
                + meshed + ": " + name);
        }
        SolverFeature f = m.sol("sol1").feature("t1");
        f.set("tout", "tsteps");
        f.set("tstepsstore", "1");
        f.set("keeplog", "on");
        f.set("maxorder", Integer.toString(maxOrder));
        f.set("rtol", Double.toString(rtol));
        f.set("tstepsbdf", "free");
        if (t0 != 0.0) f.set("starttime", num(t0));
        f.set("atolglobalmethod", method);
        f.set("atolglobalvaluemethod", valueMethod);
        if (valueMethod.equals("manual")) f.set("atolglobal", Double.toString(atolValue));
        f.set("ewtrescale", rescale);
        String[] scales = new String[variables().length];
        if (scaleMethod.equals("manual")) {
            // The variables are the Coefficient Form PDE's field first, then the
            // Global Equations variable; the specified scales follow that order.
            double[] values = model.equals("field")
                ? new double[]{scaleU, scaleV}
                : new double[]{scaleV};
            for (int i = 0; i < variables().length; ++i) {
                SolverFeature scale = m.sol("sol1").feature("v1").feature(variables()[i]);
                scale.set("scalemethod", "manual");
                scale.set("scaleval", Double.toString(values[i]));
            }
        }
        for (int i = 0; i < variables().length; ++i) {
            scales[i] = readScale(m, variables()[i]);
        }
        // The field is comp1_u and a Global Equations variable is comp1_ODE1; a
        // model without a field has no field scale to report.
        String scaleField = model.equals("field") ? scales[0] : "";
        String scaleGlobal = scales[scales.length - 1];
        f.feature("aDef").set("logsampling", "0");
        m.sol("sol1").runAll();
        double[] t = m.sol("sol1").getPVals();
        PrintWriter data = new PrintWriter(new File(directory, name + ".csv"));
        StringBuilder header = new StringBuilder("time");
        int[] size = new int[t.length];
        for (int i = 0; i < t.length; ++i) {
            size[i] = m.sol("sol1").getU(i + 1).length;
            if (i > 0 && size[i] != size[0]) {
                throw new IllegalStateException("dof count moved: " + name);
            }
        }
        for (int k = 0; k < size[0]; ++k) header.append(",u").append(k);
        for (int k = 0; k < size[0]; ++k) header.append(",d").append(k);
        data.println(header);
        for (int i = 0; i < t.length; ++i) {
            double[] u = m.sol("sol1").getU(i + 1);
            double[] d = m.sol("sol1").getUDot(i + 1);
            StringBuilder row = new StringBuilder(Double.toString(t[i]));
            for (double value : u) row.append(",").append(value);
            for (double value : d) row.append(",").append(value);
            data.println(row);
        }
        data.close();
        PrintWriter logFile = new PrintWriter(new File(directory, name + ".log"));
        logFile.print(f.getString("message"));
        logFile.close();
        String atolArray = f.getString("atol");
        String[] row = {
            name,
            model,
            num(t0),
            num(end),
            method,
            valueMethod,
            Double.toString(atolValue),
            f.getString("atolglobal"),
            f.getString("atolglobalfactor"),
            rescale,
            scaleMethod,
            Integer.toString(maxOrder),
            Integer.toString(elements),
            Double.toString(rtol),
            csv(atolArray),
            csv(f.getString("atolmethod")),
            csv(f.getString("atolvaluemethod")),
            f.getString("tderglobalmethod"),
            f.getString("tderglobalfactor"),
            csv(f.getString("tdermethod")),
            csv(f.getString("tderfactor")),
            csv(scaleField),
            csv(scaleGlobal),
            f.getString("atolglobalmethod"),
            f.getString("atolglobalvaluemethod"),
            f.getString("ewtrescale"),
            f.getString("rtol"),
            f.getString("maxorder"),
            f.getString("starttime"),
            f.getString("tout"),
            Integer.toString(size[0]),
            Double.toString(scaleU),
            Double.toString(scaleV),
        };
        manifest.println(String.join(",", row));
        manifest.flush();
        ModelUtil.clear();
        System.out.println("CASE_OK " + name + " dofs=" + size[0] + " samples=" + t.length
            + " first=" + t[0] + " second=" + t[1] + " last=" + t[t.length - 1]);
    }

    // Prints the tolerance properties once per sweep, so the names these cases
    // set can be checked against the solver instead of trusted.
    static void reportProperties() {
        Model m = build();
        SolverFeature f = m.sol("sol1").feature("t1");
        String[] properties = f.properties();
        Arrays.sort(properties);
        for (String property : properties) {
            if (!property.matches("(atol|ewt|tder|rtol|maxorder|minorder|tout|starttime).*")) {
                continue;
            }
            String value;
            try {
                value = f.getString(property);
            } catch (Exception e) {
                value = "<not a string>";
            }
            System.out.println("PROP " + property + " = " + value);
        }
        ModelUtil.clear();
    }

    public static void main(String[] args) throws Exception {
        directory = args != null && args.length > 1 ? args[1] : "tolerance";
        new File(directory).mkdirs();
        reset();
        reportProperties();
        manifest = new PrintWriter(new File(directory, "manifest.csv"));
        manifest.println("name,model,t0,end,atol_method,atol_value_method,atol_value,"
            + "atol_readback,atol_factor_readback,rescale,scale_method,max_order,elements,rtol,"
            + "atol_array,atol_method_array,atol_value_method_array,tder_method,tder_factor,"
            + "tder_method_array,tder_factor_array,scale_field,scale_global,"
            + "atol_method_readback,atol_value_method_readback,rescale_readback,"
            + "rtol_readback,max_order_readback,start_time_readback,tout_readback,export_dofs,"
            + "scale_u,scale_v");
        for (String spec : SPECS) {
            apply(spec);
            run();
        }
        manifest.close();
        System.out.println("ALL_CASES_OK");
    }
}
