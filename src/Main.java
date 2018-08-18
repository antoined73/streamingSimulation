import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class Main {

    public static void main(String[] args) {

        StreamingSimulation ss = new StreamingSimulation();
        //ss.startStreaming();

        int k_lookahead = 3;
        int deltaDl = 1;
        Double[] ps= {0.290863346226055,
                0.500000000000000,
                0.,
                0.209136653773945};
        List<Double> p_ij = Arrays.asList(ps);

        List<List<List<Double>>> s_ijl = Matrix.create3DMatrix(20,4,2,300000);
        double Ct = 1000000;
        List<List<List<Double>>> buf_it = Matrix.create3DMatrix(7, 4,2,23);
        for (List<Double> line :
                buf_it.get(0)) {
            line.set(0,1.);
        }
        Double[] jti= {2.,2.,2.,2.};
        List<Double> j_ti = Arrays.asList(jti);
        int Bmin = 2;
        int Bmax = 5;

        List<List<List<Double>>> testMat = ss.instant_optim(k_lookahead,deltaDl,p_ij,s_ijl,Ct,buf_it,j_ti,Bmin,Bmax);
        Matrix.printMatrix3D(testMat);
    }
}
