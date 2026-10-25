import com.comsol.model.*;
import com.comsol.model.util.*;
import java.io.*;
import java.util.*;

// Does COMSOL's error test average over dependent variables or over physics
// interfaces, and can a multi-variable trajectory be replayed exactly?
//
// One 1D model, built from three kinds of dependent variable:
//
//   c   A Coefficient Form PDE on the unit interval with five linear elements:
//       da=1, a=1, c=0, f=0, u(x,0) = 1+x  ->  u_i(t) = (1+x_i)*exp(-t)
//       Six solved dofs, and two internal reaction-force dofs the solver holds
//       apart from them.
//   ge  One or two Global Equations, each one dof:
//       vt + 0.5*v = 0, v(0) = 0.01     and    wt + 0.1*w = 0, w(0) = 0.1
//
// The two global equations are the discriminator: in `field_global2_one` they
// live in one Global Equations interface, in `field_global2_two` in two. The
// dependent variables, their size and their equations are identical in the two
// cases, so a rule that counts dependent variables gives the same trajectory
// and a rule that counts interfaces does not.
//
// Every case pins the variable scale to one and takes the unscaled tolerance
// criterion, so the weight of a dof is exactly `atol + rtol*abs(U_i)` and the
// record can be replayed. COMSOL's Fully Coupled node reuses its Jacobian, which
// makes its linear solve inexact and a replay impossible; each case therefore
// asks for a Jacobian rebuild at every iteration.
public class GeDofSweep {
    static final String[] SPECS = {
        // The two-variable controls: two global equations in one physics
        // interface, then the same two in an interface each.
        "field_global2_one globals=2",
        "field_global2_two globals=2 interfaces=2",
        // The multi-dof field alone: unscaled, scaled with a given scale, and
        // scaled by the solver's own.
        "field_only",
        "field_only_scaled method=scaled",
        "field_only_auto method=scaled scale=auto",
        // The closest-step storage on an output schedule that is not uniform: the
        // outputs at 0.5, 0.53 and 0.56 share one step of the run.
        "field_closest_uneven tout=closest tlist=0 0.013 0.5 0.53 0.56 2 2.05 7 10",
    };

    static String directory;
    static PrintWriter manifest;

    static String name;
    static boolean field;
    static int globals;
    static int interfaces;
    static String scaleMethod;
    static String method;
    static String tout;
    static String tlist;
    static double t0;
    static double end;
    static double rtol;
    static int maxOrder;
    static int elements;
    static boolean exact;

    static void reset() {
        name = "case";
        field = true;
        globals = 0;
        interfaces = 1;
        scaleMethod = "manual";
        method = "unscaled";
        tout = "tsteps";
        tlist = "";
        t0 = 0.0;
        end = 10.0;
        rtol = 0.01;
        maxOrder = 5;
        elements = 5;
        exact = true;
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
            if (key.equals("field")) {
                field = value.equals("on");
            } else if (key.equals("globals")) {
                globals = Integer.parseInt(value);
            } else if (key.equals("interfaces")) {
                interfaces = Integer.parseInt(value);
            } else if (key.equals("scale")) {
                scaleMethod = value;
            } else if (key.equals("method")) {
                method = value;
            } else if (key.equals("tout")) {
                tout = value;
            } else if (key.equals("tlist")) {
                // A schedule is a list, so this key takes the rest of the line.
                StringBuilder schedule = new StringBuilder(value);
                for (int j = i + 1; j < tokens.length; ++j) {
                    schedule.append(' ').append(tokens[j]);
                }
                tlist = schedule.toString();
                i = tokens.length;
            } else if (key.equals("t0")) {
                t0 = Double.parseDouble(value);
            } else if (key.equals("end")) {
                end = Double.parseDouble(value);
            } else if (key.equals("rtol")) {
                rtol = Double.parseDouble(value);
            } else if (key.equals("order")) {
                maxOrder = Integer.parseInt(value);
            } else if (key.equals("elements")) {
                elements = Integer.parseInt(value);
            } else if (key.equals("exact")) {
                exact = value.equals("on");
            } else {
                throw new IllegalArgumentException("unknown key: " + key);
            }
        }
        if (globals < 0 || globals > 2) {
            throw new IllegalArgumentException("globals must be 0, 1 or 2: " + name);
        }
        if (globals < 2 && interfaces != 1) {
            throw new IllegalArgumentException("two interfaces need two equations: " + name);
        }
        if (!scaleMethod.equals("auto") && !scaleMethod.equals("manual")) {
            throw new IllegalArgumentException("unknown scale method: " + name);
        }
        if (!method.equals("scaled") && !method.equals("unscaled")) {
            throw new IllegalArgumentException("unknown tolerance method: " + name);
        }
        if (!tout.equals("tsteps") && !tout.equals("closest")) {
            throw new IllegalArgumentException("unknown storage mode: " + name);
        }
        if (tout.equals("closest") != !tlist.isEmpty()) {
            throw new IllegalArgumentException(
                "a closest-step case needs its own schedule, and only one: " + name);
        }
        if (!field && globals == 0) {
            throw new IllegalArgumentException("a case needs a variable: " + name);
        }
    }

    static String num(double value) {
        if (value == Math.rint(value)) return Long.toString((long) value);
        return Double.toString(value);
    }

    static String csv(String value) {
        return "\"" + value.replace("\"", "\"\"") + "\"";
    }

    // The two global equations, by index: the equation, the variable name and
    // its initial value.
    static String equation(int index) {
        return index == 0 ? "vt+0.5*v" : "wt+0.1*w";
    }

    static String variableName(int index) {
        return index == 0 ? "v" : "w";
    }

    static String initialValue(int index) {
        return index == 0 ? "0.01" : "0.1";
    }

    static Model build() {
        double step = (end - t0) / 10.0;
        // The study's schedule: ten equal intervals of the span, unless the case
        // names its own, which is the output schedule of a closest-step case.
        String times = tlist.isEmpty()
            ? "range(" + num(t0) + "," + num(step) + "," + num(end) + ")"
            : tlist;
        Model m = ModelUtil.create("Model");
        m.component().create("comp1", true);
        m.component("comp1").geom().create("geom1", 1);
        m.component("comp1").geom("geom1").create("i1", "Interval");
        m.component("comp1").geom("geom1").run();
        m.component("comp1").mesh().create("mesh1");
        m.component("comp1").mesh("mesh1").feature("size").set("custom", "on");
        m.component("comp1").mesh("mesh1").feature("size").set("hmax", Double.toString(1.0 / elements));
        m.component("comp1").mesh("mesh1").run();
        if (field) {
            m.component("comp1").physics().create("c", "CoefficientFormPDE", "geom1");
            m.component("comp1").physics("c").feature("cfeq1").set("da", "1");
            m.component("comp1").physics("c").feature("cfeq1").set("a", "1");
            m.component("comp1").physics("c").feature("cfeq1").set("c", "0");
            m.component("comp1").physics("c").feature("cfeq1").set("f", "0");
            m.component("comp1").physics("c").feature("init1").set("u", "1+x");
            m.component("comp1").physics("c").prop("ShapeProperty").set("order", "1");
        }
        // One Global Equations interface holds every equation, or each equation
        // gets an interface of its own: the only difference between the two
        // cases of the discriminator. The first equation of an interface is its
        // default feature `ge1`; a second one is added to it.
        for (int i = 0; i < globals; ++i) {
            String tag = interfaces == 1 || i == 0 ? "ge" : "ge" + (i + 1);
            String feature = interfaces == 1 ? "ge" + (i + 1) : "ge1";
            if (i == 0 || interfaces != 1) {
                m.component("comp1").physics().create(tag, "GlobalEquations", "geom1");
            }
            if (!feature.equals("ge1")) {
                // The equation feature of a Global Equations interface carries the
                // interface's own type id, not a shorter name of its own.
                m.component("comp1").physics(tag).feature().create(feature, "GlobalEquations");
            }
            m.component("comp1").physics(tag).feature(feature).set("equation", equation(i));
            m.component("comp1").physics(tag).feature(feature).set("name", variableName(i));
            m.component("comp1").physics(tag).feature(feature)
                .set("initialValueU", initialValue(i));
        }
        m.study().create("std1");
        m.study("std1").create("time", "Transient");
        m.study("std1").feature("time").set("tlist", times);
        m.study("std1").createAutoSequences("all");
        return m;
    }

    // The solver's scale variables, by tag: the field first, then one per global
    // equation, whatever the solver calls them.
    static String scaleReport(Model m) {
        String[] tags = m.sol("sol1").feature("v1").feature().tags();
        StringBuilder report = new StringBuilder();
        for (String tag : tags) {
            SolverFeature feature = m.sol("sol1").feature("v1").feature(tag);
            report.append(tag).append(':').append(feature.getString("scalemethod"))
                .append(':').append(feature.getString("scaleval")).append(' ');
        }
        return report.toString().trim();
    }

    static void run() throws Exception {
        Model m = build();
        int meshed = m.component("comp1").mesh("mesh1").getNumElem();
        if (field && meshed != elements) {
            throw new IllegalStateException("expected " + elements + " elements, meshed "
                + meshed + ": " + name);
        }
        SolverFeature f = m.sol("sol1").feature("t1");
        f.set("tout", tout.equals("closest") ? "tstepsclosest" : "tsteps");
        f.set("tstepsstore", "1");
        f.set("keeplog", "on");
        f.set("maxorder", Integer.toString(maxOrder));
        f.set("rtol", Double.toString(rtol));
        f.set("tstepsbdf", "free");
        f.set("atolglobalmethod", method);
        f.feature("aDef").set("logsampling", "0");
        if (t0 != 0.0) f.set("starttime", num(t0));
        if (scaleMethod.equals("manual")) {
            for (String tag : m.sol("sol1").feature("v1").feature().tags()) {
                SolverFeature scale = m.sol("sol1").feature("v1").feature(tag);
                scale.set("scalemethod", "manual");
                scale.set("scaleval", "1");
            }
        }
        if (exact) {
            f.feature("fc1").set("jtech", "onevery");
            f.feature("fc1").set("ntolfact", "1e-8");
        }
        String jtech = f.feature("fc1").getString("jtech");
        String scales = scaleReport(m);
        m.sol("sol1").runAll();
        double[] t = m.sol("sol1").getPVals();
        int dofs = m.sol("sol1").getU(1).length;
        PrintWriter data = new PrintWriter(new File(directory, name + ".csv"));
        StringBuilder header = new StringBuilder("time");
        for (int k = 0; k < dofs; ++k) header.append(",u").append(k);
        for (int k = 0; k < dofs; ++k) header.append(",d").append(k);
        data.println(header);
        for (int i = 0; i < t.length; ++i) {
            double[] u = m.sol("sol1").getU(i + 1);
            double[] d = m.sol("sol1").getUDot(i + 1);
            if (u.length != dofs || d.length != dofs) {
                throw new IllegalStateException("dof count moved: " + name);
            }
            StringBuilder row = new StringBuilder(Double.toString(t[i]));
            for (double value : u) row.append(",").append(value);
            for (double value : d) row.append(",").append(value);
            data.println(row);
        }
        data.close();
        PrintWriter logFile = new PrintWriter(new File(directory, name + ".log"));
        logFile.print(f.getString("message"));
        logFile.close();
        String[] row = {
            name,
            field ? "on" : "off",
            Integer.toString(globals),
            Integer.toString(interfaces),
            scaleMethod,
            method,
            csv(tlist),
            num(t0),
            num(end),
            Double.toString(rtol),
            Integer.toString(maxOrder),
            Integer.toString(elements),
            exact ? "on" : "off",
            csv(jtech),
            csv(scales),
            Integer.toString(dofs),
            tout,
            f.getString("tout"),
            f.getString("atolglobalmethod"),
            f.getString("rtol"),
            f.getString("maxorder"),
            f.getString("atolglobalfactor"),
            f.getString("ewtrescale"),
            Integer.toString(meshed),
        };
        manifest.println(String.join(",", row));
        manifest.flush();
        ModelUtil.clear();
        System.out.println("CASE_OK " + name + " dofs=" + dofs + " samples=" + t.length
            + " first=" + t[0] + " second=" + t[1] + " last=" + t[t.length - 1]);
    }

    // Prints the solver properties these cases set, so the names can be checked
    // against the solver instead of trusted.
    static void reportProperties() {
        Model m = build();
        SolverFeature f = m.sol("sol1").feature("t1");
        String[] properties = f.properties();
        Arrays.sort(properties);
        for (String property : properties) {
            if (!property.matches("(atol|ewt|tder|rtol|maxorder|minorder|tout|starttime"
                + "|tsteps).*")) continue;
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
        directory = args != null && args.length > 1 ? args[1] : "dofs";
        new File(directory).mkdirs();
        reset();
        reportProperties();
        manifest = new PrintWriter(new File(directory, "manifest.csv"));
        manifest.println("name,field,globals,interfaces,scale_method,method,tlist,t0,end,rtol,"
            + "max_order,elements,exact,jtech_readback,scales_readback,dofs,tout,tout_readback,"
            + "method_readback,rtol_readback,max_order_readback,atol_factor_readback,"
            + "rescale_readback,elements_readback");
        for (String spec : SPECS) {
            apply(spec);
            run();
        }
        manifest.close();
        System.out.println("ALL_CASES_OK");
    }
}
