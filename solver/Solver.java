package solver;

import java.lang.Math;
import java.util.ArrayList;

public class Solver {

    // Coefficients (Dormand–Prince 5(4))
    final double c2 = 1.0/5.0;
    final double c3 = 3.0/10.0;
    final double c4 = 4.0/5.0;
    final double c5 = 8.0/9.0;
    
    final double a21 = 1.0/5.0;
    final double a31 = 3.0/40.0, a32 = 9.0/40.0;
    final double a41 = 44.0/45.0, a42 = -56.0/15.0, a43 = 32.0/9.0;
    final double a51 = 19372.0/6561.0, a52 = -25360.0/2187.0, a53 = 64448.0/6561.0, a54 = -212.0/729.0;
    final double a61 = 9017.0/3168.0, a62 = -355.0/33.0, a63 = 46732.0/5247.0, a64 = 49.0/176.0, a65 = -5103.0/18656.0;
    final double a71 = 35.0/384.0, a72 = 0.0, a73 = 500.0/1113.0, a74 = 125.0/192.0, a75 = -2187.0/6784.0, a76 = 11.0/84.0;

    final double b1 = a71, b2 = a72, b3 = a73, b4 = a74, b5 = a75, b6 = a76;
    
    final double bs1 = 5179.0/57600.0, bs2 = 0.0, bs3 = 7571.0/16695.0, bs4 = 393.0/640.0, 
                 bs5 = -92097.0/339200.0, bs6 = 187.0/2100.0, bs7 = 1.0/40.0;
    
    private int dimensions;
    private double q = 4.0, fac = 0.9, facmax = 5.0, facmin = 0.2;
    
    private double atol;
    private double rtol;
    
    public double t, tf, h;

    public ArrayList<double[]> y_out;
    public ArrayList<Double> t_out;
    
    private double[] k1, k2, k3, k4, k5, k6, k7, y, y4, y5, ee, sci;
    
    public Solver(double[] y0, double t0, double tf, int dimensions, double atol, double rtol) {
        this.dimensions = dimensions;
        
        // Initialize all arrays separately
        k1 = new double[dimensions];
        k2 = new double[dimensions];
        k3 = new double[dimensions];
        k4 = new double[dimensions];
        k5 = new double[dimensions];
        k6 = new double[dimensions];
        k7 = new double[dimensions];
        y4 = new double[dimensions];
        y5 = new double[dimensions];
        ee = new double[dimensions];
        sci = new double[dimensions];

        y_out = new ArrayList<>();
        t_out = new ArrayList<>();

        y_out.add(y0.clone());
        t_out.add(t0);
        
        this.y = y0.clone(); // Copy the initial values
        this.t = t0;
        this.tf = tf;
        this.atol = atol;
        this.rtol = rtol;
        this.h = initialize_stepsize(t0, y0.clone());
    }
    
    public void solve() {

        int step = 0;
        while (t < tf) {
            h = Math.min(tf - t, h);
            
            double[] y_temp = new double[dimensions];
            
            // 1st stage
            k1 = F(t, y);
            
            // 2nd stage
            for (int i = 0; i < dimensions; i++) {
                y_temp[i] = y[i] + h * (a21 * k1[i]);
            }
            k2 = F(t + c2 * h, y_temp);
            
            // 3rd stage
            for (int i = 0; i < dimensions; i++) {
                y_temp[i] = y[i] + h * (a31 * k1[i] + a32 * k2[i]);
            }
            k3 = F(t + c3 * h, y_temp);
            
            // 4th stage
            for (int i = 0; i < dimensions; i++) {
                y_temp[i] = y[i] + h * (a41 * k1[i] + a42 * k2[i] + a43 * k3[i]);
            }
            k4 = F(t + c4 * h, y_temp);
            
            // 5th stage
            for (int i = 0; i < dimensions; i++) {
                y_temp[i] = y[i] + h * (a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]);
            }
            k5 = F(t + c5 * h, y_temp);
            
            // 6th stage
            for (int i = 0; i < dimensions; i++) {
                y_temp[i] = y[i] + h * (a61 * k1[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]);
            }
            k6 = F(t + h, y_temp);
            
            // 7th stage (for embedded method)
            for (int i = 0; i < dimensions; i++) {
                y_temp[i] = y[i] + h * (a71 * k1[i] + a72 * k2[i] + a73 * k3[i] + a74 * k4[i] + a75 * k5[i] + a76 * k6[i]);
            }
            k7 = F(t + h, y_temp);
            
            // 5th order solution
            for (int i = 0; i < dimensions; i++) {
                y5[i] = y[i] + h * (b1 * k1[i] + b2 * k2[i] + b3 * k3[i] + b4 * k4[i] + b5 * k5[i] + b6 * k6[i]);
            }
            
            // 4th order embedded solution
            for (int i = 0; i < dimensions; i++) {
                y4[i] = y[i] + h * (bs1 * k1[i] + bs2 * k2[i] + bs3 * k3[i] + bs4 * k4[i] + bs5 * k5[i] + bs6 * k6[i] + bs7 * k7[i]);
            }
            
            // Error estimate
            for (int i = 0; i < dimensions; i++) {
                ee[i] = Math.abs(y5[i] - y4[i]); // Remove extra h multiplication
                sci[i] = atol + Math.max(Math.abs(y[i]), Math.abs(y5[i])) * rtol;
            }
            
            double err = hairerNorm(ee, sci);
            
            if (err <= 1.0) {
                // Step accepted
                t += h;
                y = y5.clone(); // Update solution
                y_out.add(y5.clone());
                t_out.add(t);
            }
            
            // Adjust step size for next step
            double hnew = h * Math.min(facmax, Math.max(facmin, fac * Math.pow(1.0 / Math.max(err, 1e-10), 1.0 / (q + 1.0))));
            h = hnew;
            
            step++;
            if (step > 1000) { // Safety break
                System.out.println("Too many steps, breaking");
                break;
            }
        }
		System.out.println("Steps!: " + step);
    }
    
    private double[] F(double t, double[] y) {
        double[] dydt = new double[dimensions];
        dydt[0] = y[0] - t*t + 1.0;
        return dydt;
    }

    private double hairerNorm(double[] a, double[] sci)
    {
        double _s = 0.0;
        for (int i = 0; i < dimensions; i++) {
            _s += Math.pow(a[i] / sci[i], 2.0);
        }
        return Math.sqrt(_s / dimensions);
    }

    private double initialize_stepsize(double t0, double[] y0)
    {
        // (a) Initial function evaluation and scaling
        double[] f1 = new double[dimensions];
        f1 = F(t0, y0);

        double[] sci = new double[dimensions];
        for (int i = 0; i < dimensions; i++) {
            sci[i] = atol + Math.abs(y0[i]) * rtol; // Assume m_absTol is a vector. If scalar, adjust.
        }

        double d0 = hairerNorm(y0, sci);
        double d1 = hairerNorm(f1, sci);

        // (b) First guess for the step size
        double h0;
        if (d0 < 1e-12 || d1 < 1e-12) { // Use a tighter tolerance for "zero"
            h0 = 1e-6;
        } else {
            h0 = 0.01 * (d0 / d1);
        }

        // (c) Explicit Euler step and second evaluation
        double[] y1 = new double[dimensions];
        for (int i = 0; i < dimensions; i++) {
            y1[i] = y0[i] + h0 * f1[i];
        }

        double[] f2 = new double[dimensions];
        f2 = F(t0 + h0, y1);

        // (d) Estimate the second derivative
        double[] diff_f2f1 = new double[dimensions];
        for (int i = 0; i < dimensions; i++) {
            diff_f2f1[i] = f2[i] - f1[i]; // No need for fabs here, the norm will handle it
        }
        double d2 = hairerNorm(diff_f2f1, sci) / h0;

        // (e) Refine the guess based on the error estimate
        double max_d1d2 = Math.max(d1, d2);
        double h1;

        if (max_d1d2 <= 1e-15) {
            // Handle the case where derivatives are essentially zero
            h1 = Math.max(1e-6, h0 * 1e-3);
        } else {
            double log10_h1 = ( -2.0 - Math.log10(max_d1d2) ) / (q + 1.0);
            h1 = Math.pow(10.0, log10_h1);
        }

        // (f) Final, conservative proposal
        double h = Math.min(10.0 * h0, h1);
        
        return h;
    }
}