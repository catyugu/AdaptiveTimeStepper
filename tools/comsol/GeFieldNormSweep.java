import com.comsol.model.*;
import com.comsol.model.util.*;
import java.io.*;

// Does COMSOL's time-dependent error test average per field or flat over every
// degree of freedom? One model, two independent, analytically decaying fields:
//
//   A: Global Equations, 1 dof,  uA' + rateA*uA = 0, uA(0) = 1
//   B: 1D Coefficient Form PDE on the default unit interval, linear elements,
//      da = 1, a = rateB, c = 0, f = 0, u(0) = 1   ->  B is `elements + 1`
//      dofs that all obey u' = -rateB*u exactly (M*u' + rateB*M*u = 0 with the
//      same M on both sides and a zero spatial operator), so every dof of B
//      carries the same error-to-weight ratio.
//
// With both fields the two candidate acceptance criteria differ by construction:
//   per field:  E^2 = (E_A^2 + E_B^2) / 2
//   flat:       E^2 = (E_A^2 + N_B*E_B^2) / (1 + N_B)
// where E_A and E_B are each field's own weighted RMS. Nothing else in the model
// differs between two cases of the table, so a case that changes only `elements`
// isolates the N_B factor.
//
// The order is pinned to 1 in the discriminating cases so that the ratio
// E_B/E_A = (rateB/rateA)^(p+1) is a single number for every step instead of
// following the order ramp.
//
// getU layout, measured on 6.2.0.290 rather than guessed (probes in .cache):
//   - the Global Equations dof is u[0] when the interface exists;
//   - every remaining entry is a dof of field B, except the two internal
//     reaction-force slots, which have getUDot exactly 0;
//   - the field is spatially uniform to roundoff, so one dof is the whole field
//     and `spread` (its max minus min) audits that.
//
// Every case writes name.csv (time,a,b,spread), name.log (solver log) and one
// manifest.csv row with the values read back from the solver.
public class GeFieldNormSweep {
    static final String[] SPECS = {
        // Controls: a single field, where both formulas coincide because M = 1.
        // single_a is 1 dof, single_b100 is 101 dofs (one field, uniform error).
        "single_a rateA=1",
        "single_b100 rateB=0.2236",
        "single_b25 rateB=0.2236 elements=25",
        // The sharp experiment. rateB/rateA = 0.2236 gives
        // E_B/E_A = (rateB/rateA)^2 = 0.05 at order 1.
        "two_r100 rateA=1 rateB=0.2236",
        // Same case, twice: a determinism control for the discriminating pair.
        "two_r100_repeat rateA=1 rateB=0.2236",
        // Same case at a quarter of the dofs: flat averaging must move, per-field
        // averaging must not.
        "two_r25 rateA=1 rateB=0.2236 elements=25",
        // The same mesh pair with a field that is ten times slower than A. At order
        // 1 its per-dof error ratio is (rateB/rateA)^2 = 0.01, so 101 of them still
        // weigh little against A in a flat sum: this is the pair where a flat norm
        // is diluted by the dof count and a per-field norm is not.
        "two_slow100 rateA=1 rateB=0.1",
        "two_slow25 rateA=1 rateB=0.1 elements=25",
        // Control: rateB = rateA makes the two fields identical, so both formulas
        // predict exactly the single_a sequence.
        "two_equal100 rateA=1 rateB=1",
        // The same experiment at the default order ramp, where E_B/E_A shrinks
        // as the order rises.
        "two_r100_order5 rateA=1 rateB=0.2236 order=5",
        // A second tolerance, so the conclusion does not ride on rtol = 0.01.
        "two_r100_rtol1e3 rateA=1 rateB=0.2236 rtol=1e-3",
    };

    static String directory;
    static PrintWriter manifest;

    static String name;
    static int elements;
    static double rateA;
    static double rateB;
    static int order;
    static double rtol;
    static double end;
    static boolean fullLog;

    static boolean hasGlobal;
    static boolean hasField;

    static void reset() {
        name = "case";
        elements = 100;
        rateA = 0.0;
        rateB = 0.0;
        order = 1;
        rtol = 0.01;
        end = 10.0;
        fullLog = true;
    }

    // "name key=value ..."; an omitted key keeps the default, and a rate of 0
    // omits that field from the model, so a row names only what it changes.
    static void apply(String spec) {
        reset();
        String[] tokens = spec.trim().split("\\s+");
        name = tokens[0];
        for (int i = 1; i < tokens.length; ++i) {
            int split = tokens[i].indexOf('=');
            if (split < 0) throw new IllegalArgumentException("bad token: " + tokens[i]);
            String key = tokens[i].substring(0, split);
            String value = tokens[i].substring(split + 1);
            if (key.equals("elements")) {
                elements = Integer.parseInt(value);
            } else if (key.equals("rateA")) {
                rateA = Double.parseDouble(value);
            } else if (key.equals("rateB")) {
                rateB = Double.parseDouble(value);
            } else if (key.equals("order")) {
                order = Integer.parseInt(value);
            } else if (key.equals("rtol")) {
                rtol = Double.parseDouble(value);
            } else if (key.equals("end")) {
                end = Double.parseDouble(value);
            } else if (key.equals("log")) {
                fullLog = value.equals("full");
            } else {
                throw new IllegalArgumentException("unknown key: " + key);
            }
        }
        hasGlobal = rateA > 0.0;
        hasField = rateB > 0.0;
        if (!hasGlobal && !hasField) {
            throw new IllegalArgumentException("a case needs at least one field: " + name);
        }
        if (elements < 4) throw new IllegalArgumentException("too few elements: " + name);
        if (order < 1 || order > 5) {
            throw new IllegalArgumentException("order must be 1-5: " + name);
        }
    }

    static String num(double value) {
        if (value == Math.rint(value)) return Long.toString((long) value);
        return Double.toString(value);
    }

    static Model build() {
        double step = end / 10.0;
        String times = "range(0," + num(step) + "," + num(end) + ")";
        Model m = ModelUtil.create("Model");
        m.component().create("comp1", true);
        m.component("comp1").geom().create("geom1", 1);
        // The default Interval is [0,1]; hmax = 1/elements then gives exactly
        // `elements` elements, which the run verifies from getNumElem().
        m.component("comp1").geom("geom1").create("i1", "Interval");
        m.component("comp1").geom("geom1").run();
        m.component("comp1").mesh().create("mesh1");
        m.component("comp1").mesh("mesh1").feature("size").set("custom", "on");
        m.component("comp1").mesh("mesh1").feature("size").set("hmax", Double.toString(1.0 / elements));
        m.component("comp1").mesh("mesh1").run();
        if (hasGlobal) {
            m.component("comp1").physics().create("ge", "GlobalEquations", "geom1");
            m.component("comp1").physics("ge").feature("ge1")
                .set("equation", "uAt+" + num(rateA) + "*uA");
            m.component("comp1").physics("ge").feature("ge1").set("name", "uA");
            m.component("comp1").physics("ge").feature("ge1").set("initialValueU", "1");
        }
        if (hasField) {
            m.component("comp1").physics().create("c", "CoefficientFormPDE", "geom1");
            // da*u_t + div(-c*grad(u)) + a*u = f, so c = 0 and f = 0 leave
            // u_t + rateB*u = 0 on every dof.
            m.component("comp1").physics("c").feature("cfeq1").set("da", "1");
            m.component("comp1").physics("c").feature("cfeq1").set("a", num(rateB));
            m.component("comp1").physics("c").feature("cfeq1").set("c", "0");
            m.component("comp1").physics("c").feature("cfeq1").set("f", "0");
            m.component("comp1").physics("c").feature("init1").set("u", "1");
            m.component("comp1").physics("c").prop("ShapeProperty").set("order", "1");
        }
        m.study().create("std1");
        m.study("std1").create("time", "Transient");
        m.study("std1").feature("time").set("tlist", times);
        m.study("std1").createAutoSequences("all");
        return m;
    }

    // The two fields must be weighed against the same scale, or the absolute
    // tolerance would depend on the number of dofs and the mesh comparison would
    // no longer be a comparison. The solver's automatic scale is not 1 (it read
    // back 13 for 26 field dofs and 25 for 101), so it is set explicitly.
    static String manualScale(Model m, String variable) {
        SolverFeature feature = m.sol("sol1").feature("v1").feature(variable);
        feature.set("scalemethod", "manual");
        feature.set("scaleval", "1");
        return feature.getString("scalemethod") + ":" + feature.getString("scaleval");
    }

    static String csv(String value) {
        return "\"" + value.replace("\"", "\"\"") + "\"";
    }

    static void run() throws Exception {
        Model m = build();
        int meshed = m.component("comp1").mesh("mesh1").getNumElem();
        if (meshed != elements) {
            throw new IllegalStateException("expected " + elements + " elements, meshed "
                + meshed + ": " + name);
        }
        SolverFeature f = m.sol("sol1").feature("t1");
        f.set("tout", "tsteps");
        f.set("tstepsstore", "1");
        f.set("keeplog", "on");
        f.set("maxorder", Integer.toString(order));
        f.set("rtol", Double.toString(rtol));
        f.set("tstepsbdf", "free");
        String scaleU = "";
        String scaleA = "";
        if (hasField) scaleU = manualScale(m, "comp1_u");
        if (hasGlobal) scaleA = manualScale(m, "comp1_ODE1");
        if (fullLog) f.feature("aDef").set("logsampling", "0");
        String elementOrder = hasField
            ? m.component("comp1").physics("c").prop("ShapeProperty").getString("order") : "";
        m.sol("sol1").runAll();
        PrintWriter data = new PrintWriter(new File(directory, name + ".csv"));
        data.println("time,a,b,spread");
        double[] t = m.sol("sol1").getPVals();
        int fieldDofs = 0;
        for (int i = 0; i < t.length; ++i) {
            double[] u = m.sol("sol1").getU(i + 1);
            double[] d = m.sol("sol1").getUDot(i + 1);
            double a = hasGlobal ? u[0] : Double.NaN;
            int first = hasGlobal ? 1 : 0;
            double b = Double.NaN;
            double low = Double.POSITIVE_INFINITY;
            double high = Double.NEGATIVE_INFINITY;
            int seen = 0;
            for (int k = first; k < u.length; ++k) {
                // The two internal reaction-force slots never move.
                if (d[k] == 0.0) continue;
                if (seen == 0) b = u[k];
                low = Math.min(low, u[k]);
                high = Math.max(high, u[k]);
                ++seen;
            }
            if (hasField) {
                if (fieldDofs == 0) {
                    fieldDofs = seen;
                } else if (seen != fieldDofs) {
                    throw new IllegalStateException("field dof count moved: " + name);
                }
            }
            data.println(Double.toString(t[i]) + "," + Double.toString(a) + ","
                + Double.toString(b) + ","
                + Double.toString(hasField ? high - low : Double.NaN));
        }
        data.close();
        if (hasField && fieldDofs != elements + 1) {
            throw new IllegalStateException("expected " + (elements + 1) + " field dofs, found "
                + fieldDofs + ": " + name);
        }
        PrintWriter logFile = new PrintWriter(new File(directory, name + ".log"));
        logFile.print(f.getString("message"));
        logFile.close();
        String[] row = {
            name,
            Integer.toString(elements),
            Integer.toString(meshed),
            Integer.toString(fieldDofs),
            hasGlobal ? "global" : "none",
            hasField ? "pde" : "none",
            Double.toString(rateA),
            Double.toString(rateB),
            Integer.toString(order),
            Double.toString(rtol),
            Double.toString(end),
            fullLog ? "full" : "default",
            csv(scaleU),
            csv(scaleA),
            csv(elementOrder),
            f.getString("maxorder"),
            f.getString("rtol"),
            f.getString("tout"),
        };
        manifest.println(String.join(",", row));
        manifest.flush();
        ModelUtil.clear();
        System.out.println("CASE_OK " + name + " field_dofs=" + fieldDofs
            + " samples=" + t.length + " last=" + t[t.length - 1]);
    }

    public static void main(String[] args) throws Exception {
        directory = args != null && args.length > 1 ? args[1] : "twofield";
        new File(directory).mkdirs();
        manifest = new PrintWriter(new File(directory, "manifest.csv"));
        manifest.println("name,elements,elements_readback,field_dofs,field_a,field_b,"
            + "rate_a,rate_b,max_order,rtol,end,log_sampling,field_scale,global_scale,"
            + "element_order,max_order_readback,rtol_readback,tout_readback");
        for (String spec : SPECS) {
            apply(spec);
            run();
        }
        manifest.close();
        System.out.println("ALL_CASES_OK");
    }
}
