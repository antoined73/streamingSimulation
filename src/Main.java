import java.io.IOException;
import java.util.Arrays;
import java.util.List;

public class Main {

    public static void main(String[] args) {

        List<List<List<Double>>> s_ijl = null;
        try {

            CSVReader csvReader = new CSVReader("./sample.csv",
                    4,
                    177,
                    5);

            s_ijl = csvReader.read();
        } catch (IOException e) {
            e.printStackTrace();
        }
        System.out.println("s_ijl");
        Matrix.printMatrix3D(s_ijl);
        System.out.println("------------");
        StreamingSimulation ss = new StreamingSimulation(s_ijl);
        ss.startStreaming();

        /***
        int k_lookahead = 3;
        int deltaDl = 1;
        Double[] ps= {0.342535707281441,
                0.,
                0.500000000000000,
                0.157464292718559};
        List<Double> p_ij = Arrays.asList(ps);

        CSVWriter csvWriter = new CSVWriter();

        List<List<List<Double>>> s_ijl = null;
        try {
            csvWriter.write();

            System.out.println(csvWriter.qualities);
            CSVReader csvReader = new CSVReader("./sample.csv",
                    (1+csvWriter.maxX) * (1+csvWriter.maxY),
                    csvWriter.maxSegmentNumber,
                    csvWriter.qualities.size());

            s_ijl = csvReader.read();
        } catch (IOException e) {
            e.printStackTrace();
        }

        Matrix.printMatrix3D(s_ijl);


        List<List<List<Double>>> s_ijl = Matrix.create3DMatrix(20,4,2,150000);
        for (List<Double> v :
                s_ijl.get(1)) {
            for (int i=0; i<v.size(); i++) {
                v.set(i, v.get(i)*2);
            }
        }

        double Ct = 1000000;
        List<List<List<Double>>> buf_it = Matrix.create3DMatrix(7, 4, 2, -3);
        for (List<Double> line :
                buf_it.get(0)) {
            line.set(0,1.);
        }
        Double[] jti= {2.,2.,2.,2.};
        List<Double> j_ti = Arrays.asList(jti);
        int Bmin = 2;
        int Bmax = 5;
        double time_video = 0.;


        List<List<List<Double>>> testMat = ss.instant_optim(k_lookahead,deltaDl,p_ij,s_ijl,Ct,buf_it,j_ti,Bmin,Bmax,time_video);
        Matrix.printMatrix3D(testMat);***/
    }
}
