import com.comsol.model.*;
import com.comsol.model.util.*;
import java.io.*;

// The same second-order ODE in COMSOL's two formulations:
//
//   second_order      one Global Equation written in qtt: qtt+0.1*qt+q=0
//   first_order_pair  two Global Equations: at-b=0, bt+0.1*b+a=0
//   pair_scale1       the pair with both variable scales pinned to one
//
// q(0)=1, q'(0)=100. The derivative is two orders of magnitude above the state,
// so an error test that covers the auxiliary velocity computes a norm the state
// alone does not: the two readings differ by sqrt(2) in the first step.
//
// Every case writes name.csv (time, all dofs, all dof derivatives), name.log (the
// solver log, which carries the scale block and the accepted-step table) and one
// manifest.csv row. The record of the second-order form is
// test/reference/comsol_second_order.hpp; the two pair runs are the controls
// quoted in doc/evidence/second-order.md.
public class GeSecondOrder {
    static String directory;
    static PrintWriter manifest;

    static Model build(boolean secondOrder) {
        Model m = ModelUtil.create("Model");
        m.component().create("comp1", true);
        m.component("comp1").geom().create("geom1", 1);
        m.component("comp1").geom("geom1").create("i1", "Interval");
        m.component("comp1").geom("geom1").run();
        m.component("comp1").mesh().create("mesh1");
        m.component("comp1").mesh("mesh1").run();
        m.component("comp1").physics().create("ge", "GlobalEquations", "geom1");
        if (secondOrder) {
            m.component("comp1").physics("ge").feature("ge1").set("equation", "qtt+0.1*qt+q");
            m.component("comp1").physics("ge").feature("ge1").set("name", "q");
            m.component("comp1").physics("ge").feature("ge1").set("initialValueU", "1");
            m.component("comp1").physics("ge").feature("ge1").set("initialValueUt", "100");
        } else {
            m.component("comp1").physics("ge").feature("ge1").set("equation", "at-b");
            m.component("comp1").physics("ge").feature("ge1").set("name", "a");
            m.component("comp1").physics("ge").feature("ge1").set("initialValueU", "1");
            m.component("comp1").physics("ge").create("ge2", "GlobalEquations");
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

    static void run(String name, boolean secondOrder, boolean unitScale) throws Exception {
        Model m = build(secondOrder);
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
        // Every accepted step in the log, which is where the step sizes and the
        // Tfail counter come from.
        f.feature("aDef").set("logsampling", "0");
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
            Integer.toString(t.length), Double.toString(t[0]), Double.toString(t[1]),
            Double.toString(t[t.length - 1]), "\"" + scales.replace("\"", "\"\"") + "\""));
        manifest.flush();
        ModelUtil.clear();
        System.out.println("CASE_OK " + name + " dofs=" + dofs + " samples=" + t.length
            + " second=" + t[1] + " last=" + t[t.length - 1] + " scales=[" + scales + "]");
    }

    public static void main(String[] args) throws Exception {
        directory = args != null && args.length > 1 ? args[1] : "secondorder";
        new File(directory).mkdirs();
        manifest = new PrintWriter(new File(directory, "manifest.csv"));
        manifest.println("name,dofs,samples,first,second,last,scales");
        run("second_order", true, false);
        run("first_order_pair", false, false);
        run("pair_scale1", false, true);
        manifest.close();
        System.out.println("ALL_CASES_OK");
    }
}