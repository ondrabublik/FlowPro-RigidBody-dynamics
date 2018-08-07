package flowpro.user.dynamics;

import flowpro.user.auxiliary.ScriptEvaluator;
import flowpro.api.*;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import javax.script.ScriptException;

/**
 *
 * @author obublik
 */
public class Rigid2D implements Dynamics {

    Equation eqn;

    private int nBodies;
    private Body[] bodies;
    private String simulationPath;
    private double[] xForce;
    private double[] yForce;
    private double[] momentum;
    public double dt;
    public double t;
    public double zLength;
    public double tKick;

    // dynamic
    boolean dynamicComputation = false;
    protected double lRef = 1;
    protected double pRef = 1;
    protected double rhoRef = 1;
    protected double vRef = 1;
    protected double tRef = 1;
    protected double mRef;
    protected double IRef;
    protected double kRef;
    protected double torRef;
    protected double bRef;
    protected double torbRef;

    public void init(int nBodies, String simulationPath, String meshPath, Equation eqn) throws IOException {
        this.eqn = eqn;
        this.nBodies = nBodies;
        this.simulationPath = simulationPath;
        ScriptEvaluator jsEval = new ScriptEvaluator();
        bodies = new Body[nBodies];
        for (int i = 0; i < nBodies; i++) {
            bodies[i] = new Body(i, jsEval);
        }

        // clear file for results
        FlowProProperties props = new FlowProProperties();
        props.load(new FileInputStream(simulationPath + "parameters.txt"));
        if (props.containsKey("continueComputation")) {
            try {
                PrintWriter writer = new PrintWriter(simulationPath + "bodiesDynamic.txt");
                writer.print("");
                writer.close();
            } catch (IOException ioe) {
                System.out.println("Error while creating a new empty file :" + ioe);
            }
            try {
                PrintWriter writer = new PrintWriter(simulationPath + "bodiesUserDef.txt");
                writer.print("");
                writer.close();
            } catch (IOException ioe) {
                System.out.println("Error while creating a new empty file :" + ioe);
            }
        }

        // try to read body centers
        try {
            double[] rotationCenters = props.getDoubleArray("rotationCenters");
            for (int i = 0; i < 2 * nBodies; i = i + 2) {
                bodies[i].XCenter = new double[]{rotationCenters[i], rotationCenters[i + 1]}; //???? chyba
            }
        } catch (IOException ioe) {
            System.out.println("Body centers not defined!");
        }

        zLength = 1;
        if (props.containsKey("zLength")) {
            zLength = props.getDouble("zLength");
        } else {
            System.out.println("Length of 2D bodies in z coordinate is set to " + zLength);
        }

        tKick = Double.MAX_VALUE;
        if (props.containsKey("tKick")) {
            tKick = props.getDouble("tKick");
        }

        try {
            double[] refValues = eqn.getReferenceValues();
            lRef = refValues[0];
            pRef = refValues[1];
            rhoRef = refValues[2];
            tRef = refValues[4];
        } catch (Exception e) {
            System.out.println("Cannot assign referential values to body dynamics!");
        }

        try {
            double[] m = props.getDoubleArray("M");
            double[] b = props.getDoubleArray("B");
            double[] k = props.getDoubleArray("K");
            int n = (int) Math.sqrt(m.length);
            double[][] M = new double[n][n];
            double[][] B = new double[n][n];
            double[][] K = new double[n][n];
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    M[i][j] = m[n * i + j];
                    B[i][j] = b[n * i + j];
                    K[i][j] = k[n * i + j];
                }
            }

            /*mRef = rhoRef * lRef * lRef;
             kRef = pRef;
             bRef = pRef * tRef;
             IRef = mRef * lRef * lRef;
             torRef = pRef * lRef * lRef;
             torbRef = pRef * lRef * lRef * tRef;
             //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! referencni hodnoty deviacnich momentu!!!!!!!!!!!
             double[][] MRef = new double[][]{{1 / mRef, 0, 0}, {0, 1 / mRef, 0}, {0, 0, 1 / IRef}};
             double[][] BRef = new double[][]{{1 / bRef, 0, 0}, {0, 1 / bRef, 0}, {0, 0, 1 / torbRef}};
             double[][] KRef = new double[][]{{1 / kRef, 0, 0}, {0, 1 / kRef, 0}, {0, 0, 1 / torRef}};

             // transfer to unreferential values
             M = Mat.times(M, MRef);
             B = Mat.times(B, BRef);
             K = Mat.times(K, KRef);*/
            double[][] iFRef = new double[][]{{1 / (pRef * lRef), 0, 0}, {0, 1 / (pRef * lRef), 0}, {0, 0, 1 / (pRef * lRef * lRef)}};
            double[][] LRef = new double[][]{{lRef, 0, 0}, {0, lRef, 0}, {0, 0, lRef * lRef}};
            // transfer to unreferential values
            M = Mat.times(Mat.times(Mat.times(iFRef, M), LRef), 1 / (tRef * tRef));
            B = Mat.times(Mat.times(Mat.times(iFRef, B), LRef), 1 / tRef);
            K = Mat.times(Mat.times(iFRef, K), LRef);

            for (int i = 0; i < nBodies; i++) {
                bodies[i].setBodyParameters(M, B, K);
            }

            dynamicComputation = true;

        } catch (Exception e) {
            System.out.println("Mass, stifness or damping matrix not defined!");
        }

        // load scripts
        String[] xMotion = null;
        String[] yMotion = null;
        String[] alphaMotion = null;
        try {
            if (props.containsKey("xMotion")) {
                xMotion = props.getStringArray("xMotion");
            }
            if (props.containsKey("yMotion")) {
                yMotion = props.getStringArray("yMotion");
            }
            if (props.containsKey("alphaMotion")) {
                alphaMotion = props.getStringArray("alphaMotion");
            }
        } catch (Exception e) {
            System.out.println("Error in kinematic scripts " + e);
        }
        for (int i = 0; i < nBodies; i++) {
            if (xMotion != null) {
                bodies[i].setXMotion(xMotion[i]);
            }
            if (yMotion != null) {
                bodies[i].setYMotion(yMotion[i]);
            }
            if (alphaMotion != null) {
                bodies[i].setAlphaMotion(alphaMotion[i]);
            }
        }

        String[] xTimeForce = null;
        String[] yTimeForce = null;
        String[] alphaTimeForce = null;
        try {
            if (props.containsKey("xExternForceInTime")) {
                xTimeForce = props.getStringArray("xExternForceInTime");
            }
            if (props.containsKey("yExternForceInTime")) {
                yTimeForce = props.getStringArray("yExternForceInTime");
            }
            if (props.containsKey("alphaExternForceInTime")) {
                alphaTimeForce = props.getStringArray("alphaExternForceInTime");
            }
        } catch (Exception e) {
            System.out.println("Error in kinematic scripts " + e);
        }
        for (int i = 0; i < nBodies; i++) {
            if (xTimeForce != null) {
                bodies[i].setXTimeForce(xTimeForce[i]);
            }
            if (yTimeForce != null) {
                bodies[i].setYTimeForce(yTimeForce[i]);
            }
            if (alphaTimeForce != null) {
                bodies[i].setAlphaTimeForce(alphaTimeForce[i]);
            }
        }
    }

    public void computeBodyMove(double dt, double t, FluidForces fluFor) {
        this.dt = dt;
        this.t = t;
        double[][] translationForce = fluFor.getTranslationForce();
        this.xForce = translationForce[0];
        this.yForce = translationForce[1];
        double[][] rotationForce = fluFor.getRotationForce();
        this.momentum = rotationForce[0];

        try {
            if (dynamicComputation) {
                for (int i = 0; i < nBodies; i++) {
                    bodies[i].F = new double[]{zLength * xForce[i], zLength * yForce[i], zLength * momentum[i]};
                    bodies[i].F = Mat.plusVec(bodies[i].F, bodies[i].timeDependentForces(t));
                    double[] BU = Mat.times(bodies[i].B, bodies[i].U);
                    double[] KX = Mat.times(bodies[i].K, bodies[i].X);
                    //bodies[i].RHS = Mat.times(bodies[i].iM, Mat.plusMinMinVec(bodies[i].F, BU, KX));
                    double[] aux = new double[3];
                    for (int j = 0; j < 3; j++) {
                        aux[j] = -(BU[j] + KX[j]);
                    }
                    bodies[i].RHS = Mat.times(bodies[i].iM, aux);
                    // Two-step Adams–Bashforth
                    for (int j = 0; j < bodies[i].X.length; j++) {
                        bodies[i].Xnew[j] = bodies[i].X[j] + dt * (1.5 * bodies[i].U[j] - 0.5 * bodies[i].Uold[j]);
                        bodies[i].Unew[j] = bodies[i].U[j] + dt * (1.5 * bodies[i].RHS[j] - 0.5 * bodies[i].RHSold[j]) + dt * (bodies[i].F[j] + bodies[i].Fold[j]) / 2;
                    }
                }
            }

            // kinematic forced body
            if (t < tKick) {
                for (int i = 0; i < nBodies; i++) {
                    bodies[i].setActualKinematicCoordinates(t);
                }
            }
        } catch (Exception e) {
            System.out.println("formal error in JavaScript script ");
        }
    }

    public void nextTimeLevel() {
        if (dynamicComputation) {
            for (int i = 0; i < nBodies; i++) {
                for (int j = 0; j < bodies[i].X.length; j++) {
                    bodies[i].Uold[j] = bodies[i].U[j];
                    bodies[i].RHSold[j] = bodies[i].RHS[j];
                    bodies[i].Fold[j] = bodies[i].F[j];
                    bodies[i].X[j] = bodies[i].Xnew[j];
                    bodies[i].U[j] = bodies[i].Unew[j];
                }
            }
        }
    }

    public MeshMove[] getMeshMove() {
        MeshMove[] mshMov = new MeshMove[nBodies];
        for (int k = 0; k < nBodies; k++) {
            mshMov[k] = new MeshMove(new double[]{bodies[k].Xnew[0], bodies[k].Xnew[1]}, new double[]{bodies[k].Xnew[2]}, null, null);
        }
        return mshMov;
    }

    public double[][] getCenter() { // XCenter neni updatovan ?????????
        double[][] center = new double[2][nBodies];
        for (int i = 0; i < nBodies; i++) {
            center[0][i] = bodies[i].XCenter[0];
            center[1][i] = bodies[i].XCenter[1];
        }
        return center;
    }

    public void savePositionsAndForces() {
        try (FileWriter fw = new FileWriter(simulationPath + "bodiesDynamic.txt", true);
                BufferedWriter bw = new BufferedWriter(fw);
                PrintWriter out = new PrintWriter(bw)) {
            String line = Double.toString(t);
            for (int i = 0; i < nBodies; i++) {
                line = line + " " + Double.toString(bodies[i].Xnew[0]) + " " + Double.toString(bodies[i].Xnew[1]) + " " + Double.toString(bodies[i].Xnew[2]) + " "
                        + Double.toString(xForce[i]) + " " + Double.toString(yForce[i]) + " " + Double.toString(momentum[i]);
            }
            out.println(line);
        } catch (IOException e) {
            //exception handling left as an exercise for the reader
        }
    }

    public class Body {

        ScriptEvaluator jsEval;
        String xMotion;
        String yMotion;
        String alphaMotion;
        String xTimeForce;
        String yTimeForce;
        String alphaTimeForce;

        public double[] X;
        public double[] U;

        public double[][] M;
        public double[][] iM;
        public double[][] B;
        public double[][] K;

        public double[] XCenter;

        // old values
        public double[] Xnew;
        public double[] Unew;
        public double[] Uold;
        public double[] RHS;
        public double[] RHSold;
        public double[] F;
        public double[] Fold;

        Body(int i, ScriptEvaluator jsEval) {
            X = new double[3];
            U = new double[3];
            XCenter = new double[3];

            Xnew = new double[3];
            Unew = new double[3];
            Uold = new double[3];
            RHS = new double[3];
            RHSold = new double[3];
            Fold = new double[3];

            this.jsEval = jsEval;
        }

        void setBodyParameters(double[][] M, double[][] B, double[][] K) {
            this.M = M;
            this.B = B;
            this.K = K;
            iM = Mat.invert(M);
        }

        void setActualKinematicCoordinates(double t) throws ScriptException {
            t = t * tRef;
            if (xMotion != null) {
                Xnew[0] = jsEval.eval(xMotion, t);
            }
            if (yMotion != null) {
                Xnew[1] = jsEval.eval(yMotion, t);
            }
            if (alphaMotion != null) {
                Xnew[2] = jsEval.eval(alphaMotion, t);
            }
        }

        public double[] timeDependentForces(double t) throws ScriptException { // time dependent external forces
            t = t * tRef;
            double[] Ft = new double[3];
            if (xTimeForce != null) {
                Ft[0] = jsEval.eval(xTimeForce, t);
            }
            if (yTimeForce != null) {
                Ft[1] = jsEval.eval(yTimeForce, t);
            }
            if (alphaTimeForce != null) {
                Ft[2] = jsEval.eval(alphaTimeForce, t);
            }
            return Ft;
        }

        public void setXMotion(String xMotion) {
            this.xMotion = xMotion;
        }

        public void setYMotion(String yMotion) {
            this.yMotion = yMotion;
        }

        public void setAlphaMotion(String alphaMotion) {
            this.alphaMotion = alphaMotion;
        }

        public void setXTimeForce(String xTimeForce) {
            this.xTimeForce = xTimeForce;
        }

        public void setYTimeForce(String yTimeForce) {
            this.yTimeForce = yTimeForce;
        }

        public void setAlphaTimeForce(String alphaTimeForce) {
            this.alphaTimeForce = alphaTimeForce;
        }
    }
}
