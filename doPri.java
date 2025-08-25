
import solver.Solver;

public class doPri {
    public static void main(String[] args) {
        double[] y0 = new double[1];
        y0[0] = 0.5;
        double t0 = 0.0;
        double tf = 2.0;
        double h = 1e-3;
        double aTol = 1e-6;
        double rTol = 1e-4;

        int dim = 1;
        
        Solver solver = new Solver(y0, t0, tf, dim, aTol, rTol);
        solver.solve();

        for (int i=0; i<solver.y_out.size(); i++)
        {
            System.out.println("t: " + solver.t_out.get(i));
            for (int j=0; j<dim ; j++)
            {
                double[] sol = solver.y_out.get(i);
                System.out.print(" " + sol[j]);
            }
            System.out.println();
        }
    }
}