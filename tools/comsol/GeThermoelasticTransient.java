import com.comsol.model.Model;
import com.comsol.model.SolverFeature;
import com.comsol.model.util.ModelUtil;
import java.io.PrintWriter;

// Nondimensional reciprocal thermoelasticity, independently solved by COMSOL.
// q0=x, q1=v, q2=theta. Three eigenmodes have decay rates 1, 2 and 3.
public class GeThermoelasticTransient {
    public static void main(String[] args) throws Exception {
        Model m = ModelUtil.create("Model");
        m.component().create("comp1", true);
        m.component("comp1").geom().create("geom1", 2);
        m.component("comp1").geom("geom1").create("sq1", "Square");
        m.component("comp1").geom("geom1").run();
        m.component("comp1").physics().create("ge", "GlobalEquations", "geom1");
        String[] equations = {"q0t-q1", "q1t+4*q0+4.5*q1-0.5*q2", "q2t+0.5*q1+1.5*q2"};
        String[] initial = {"1.75", "-2.75", "-0.25"};
        for (int i = 0; i < 3; ++i) {
            String tag = "ge" + (i + 1);
            if (i > 0)
                m.component("comp1").physics("ge").create(tag, "GlobalEquations");
            m.component("comp1").physics("ge").feature(tag).set("name", "q" + i);
            m.component("comp1").physics("ge").feature(tag).set("equation", equations[i]);
            m.component("comp1").physics("ge").feature(tag).set("initialValueU", initial[i]);
        }
        m.component("comp1").mesh().create("mesh1");
        m.component("comp1").mesh("mesh1").autoMeshSize(9);
        m.component("comp1").mesh("mesh1").run();
        m.study().create("std1");
        m.study("std1").create("time", "Transient");
        m.study("std1").feature("time").set("tlist", "0 5");
        m.study("std1").createAutoSequences("all");
        SolverFeature f = m.sol("sol1").feature("t1");
        f.set("tout", "tsteps");
        f.set("keeplog", "on");
        f.set("maxorder", 5);
        f.set("rtol", 1e-7);
        f.set("atolglobalmethod", "unscaled");
        f.set("atolglobalvaluemethod", "factor");
        f.set("atolglobalfactor", 0.1);
        f.set("bwinitstepfrac", 1e-5);
        f.set("endtimeinterpolation", "off");
        f.feature("fc1").set("jtech", "onevery");
        f.feature("fc1").set("ntolfact", 1e-8);
        m.sol("sol1").runAll();

        double[] times = m.sol("sol1").getPVals();
        PrintWriter data = new PrintWriter("thermoelastic.csv");
        data.println("time,x,v,theta,dx,dv,dtheta");
        for (int i = 0; i < times.length; ++i) {
            double[] y = m.sol("sol1").getU(i + 1);
            double[] dy = m.sol("sol1").getUDot(i + 1);
            if (y.length != 3 || dy.length != 3)
                throw new IllegalStateException("expected three global degrees of freedom");
            data.print(Double.toString(times[i]));
            for (double value : y) data.print("," + Double.toString(value));
            for (double value : dy) data.print("," + Double.toString(value));
            data.println();
        }
        data.close();
        PrintWriter log = new PrintWriter("thermoelastic.log");
        log.print(f.getString("message"));
        log.close();
        m.save("thermoelastic.mph");
        System.out.println("THERMOELASTIC_OK samples=" + times.length);
    }
}
