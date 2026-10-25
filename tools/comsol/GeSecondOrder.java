import com.comsol.model.*;
import com.comsol.model.util.*;
import java.io.*;

// The same second-order ODE in COMSOL's formulations, and copies of it:
//
//   second_order       one Global Equation written in qtt: qtt+0.1*qt+q=0
//   second_order_exact the same with its Jacobian rebuilt at every iteration
//   second_order_k2    two independent copies of it, one variable each
//   second_order_k3    three independent copies
//   first_order_pair   two Global Equations: at-b=0, bt+0.1*b+a=0
//   pair_scale1        the pair with both variable scales pinned to one
//
// q(0)=1, q'(0)=100. The derivative is two orders of magnitude above the state,
// so an error test that covers the auxiliary velocity computes a norm the state
// alone does not.
//
// The copies ask how many variables the error test's average divides by: with
// K independent oscillators the per-variable normalized error is the same one
// K times, so a divisor of 2K leaves the first step alone, a divisor of K + 1
// shortens it by sqrt((K+1)/(2K)), a divisor of K by sqrt(2), and a divisor of
// two by sqrt(K).
//
// Every case writes name.csv (time, all dofs, all dof derivatives), name.log (the
// solver log, which carries the scale block and the accepted-step table) and one
// manifest.csv row. The records are test/reference/comsol_second_order.hpp and
// test/reference/comsol_second_order_k.hpp; the two pair runs are the controls
// quoted in doc/evidence/second-order.md.
public class GeSecondOrder {
    static String directory;
    static PrintWriter manifest;

    static Model build(boolean secondOrder, int oscillators) {
        Model m = ModelUtil.create("Model");
        m.component().create("comp1", true);
        m.component("comp1").geom().create("geom1", 1);
        m.component("comp1").geom("geom1").create("i1", "Interval");
        m.component("comp1").geom("geom1").run();
        m.component("comp1").mesh().create("mesh1");
        m.component("comp1").mesh("mesh1").run();
        m.component("comp1").physics().create("ge", "GlobalEquations", "geom1");
        if (secondOrder) {
            for (int k = 1; k <= oscillators; ++k) {
                String name = k == 1 ? "q" : "q" + k;
                String tag = k == 1 ? "ge1" : "ge" + k;
                if (k > 1) {
                    m.component("comp1").physics("ge").feature().create(tag, "GlobalEquations");
                }
                m.component("comp1").physics("ge").feature(tag)
                    .set("equation", name + "tt+0.1*" + name + "t+" + name);
                m.component("comp1").physics("ge").feature(tag).set("name", name);
                m.component("comp1").physics("ge").feature(tag).set("initialValueU", "1");
                m.component("comp1").physics("ge").feature(tag).set("initialValueUt", "100");
            }
        } else {
            m.component("comp1").physics("ge").feature("ge1").set("equation", "at-b");
            m.component("comp1").physics("ge").feature("ge1").set("name", "a");
            m.component("comp1").physics("ge").feature("ge1").set("initialValueU", "1");
            m.component("comp1").physics("ge").feature().create("ge2", "GlobalEquations");
            m.component("comp1").physics("ge").feature("ge2").set("equation", "bt+0.1*b+a");
            m.component("comp1").physics("ge").feature("ge2").set("name", "b");
            m.component("comp1").physics("ge").feature("ge2").set("initialValueU", "100");
        }
        m.study().create("std1");
        m.study("std1").create("time", "Transient");
        m.study("std1").feature("time").set("tlist", "range(0,1,10)");
        m.study("std1").createAutoSequences("all");
        return m;
    }

    static void run(String name, boolean secondOrder, boolean unitScale, boolean exact,
        int oscillators) throws Exception {
        Model m = build(secondOrder, oscillators);
        if (unitScale) {
            // The pair pinned to the scale the second-order form reports for its
            // own variable: one. What is left between the formulations is then
            // only which components the error test covers.
            String[] variables = {"comp1_ODE1", "comp1_ODE2"};
            for (int i = 0; i < variables.length; ++i) {
                m.sol("sol1").feature("v1").feature(variables[i]).set("scalemethod", "manual");
                m.sol("sol1").feature("v1").feature(variables[i]).set("scaleval", "1");
            }
        }
        SolverFeature f = m.sol("sol1").feature("t1");
        f.set("tout", "tsteps");
        f.set("keeplog", "on");
        f.set("maxorder", "5");
        f.set("rtol", "0.01");
        f.set("tstepsbdf", "free");
        f.set("atolglobalmethod", "scaled");
        f.set("atolglobalvaluemethod", "factor");
        if (exact) {
            // The default Fully Coupled node reuses its Jacobian, so a linear
            // model's solve stops when the increment is small, not when the
            // residual is. This case takes the exact root instead.
            f.feature("fc1").set("jtech", "onevery");
            f.feature("fc1").set("ntolfact", "1e-8");
        }
        // Every accepted step in the log, which is where the step sizes and the
        // Tfail counter come from.
        f.feature("aDef").set("logsampling", "0");
        m.sol("sol1").runAll();
        double[] t = m.sol("sol1").getPVals();
        int dofs = m.sol("sol1").getU(1).length;
        if (secondOrder && dofs != oscillators) {
            throw new IllegalStateException("expected " + oscillators + " dofs, got " + dofs);
        }
        PrintWriter data = new PrintWriter(new File(directory, name + ".csv"));
        StringBuilder header = new StringBuilder("time");
        for (int k = 0; k < dofs; ++k) header.append(",u").append(k);
        for (int k = 0; k < dofs; ++k) header.append(",d").append(k);
        data.println(header);
        for (int i = 0; i < t.length; ++i) {
            double[] u = m.sol("sol1").getU(i + 1);
            double[] d = m.sol("sol1").getUDot(i + 1);
            if (u.length != dofs) throw new IllegalStateException("dof count moved: " + name);
            StringBuilder row = new StringBuilder(Double.toString(t[i]));
            for (double value : u) row.append(",").append(value);
            for (double value : d) row.append(",").append(value);
            data.println(row);
        }
        data.close();
        String log = f.getString("message");
        PrintWriter logFile = new PrintWriter(new File(directory, name + ".log"));
        logFile.print(log);
        logFile.close();
        String scales = "";
        int block = log.indexOf("Scales for dependent variables:");
        if (block >= 0) {
            String[] lines = log.substring(block).split("\n");
            for (int i = 1; i < lines.length; ++i) {
                String line = lines[i].trim();
                if (line.length() == 0 || line.startsWith("Step")) break;
                if (scales.length() > 0) scales = scales + " ";
                scales = scales + line;
            }
        }
        manifest.println(String.join(",", name, Integer.toString(dofs),
            Integer.toString(oscillators), Integer.toString(t.length),
            Double.toString(t[0]), Double.toString(t[1]), Double.toString(t[t.length - 1]),
            exact ? "onevery" : "default", "\"" + scales.replace("\"", "\"\"") + "\""));
        manifest.flush();
        ModelUtil.clear();
        System.out.println("CASE_OK " + name + " dofs=" + dofs + " samples=" + t.length
            + " second=" + t[1] + " h1=" + (t[1] - 2e-5) + " last=" + t[t.length - 1]
            + " scales=[" + scales + "]");
    }

    public static void main(String[] args) throws Exception {
        directory = args != null && args.length > 1 ? args[1] : "secondorder";
        new File(directory).mkdirs();
        manifest = new PrintWriter(new File(directory, "manifest.csv"));
        manifest.println("name,dofs,oscillators,samples,first,second,last,jtech,scales");
        run("second_order", true, false, false, 1);
        run("second_order_exact", true, false, true, 1);
        run("second_order_k2", true, false, true, 2);
        run("second_order_k3", true, false, true, 3);
        run("first_order_pair", false, false, false, 0);
        run("pair_scale1", false, true, false, 0);
        manifest.close();
        System.out.println("ALL_CASES_OK");
    }
}
