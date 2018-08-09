import java.awt.*;
import java.util.ArrayList;
import java.util.List;

public class Matrix {


    public static void multiply3D(double multiplicator, List<List<List<Double>>> matrix3D) {
        for(int i=0; i<matrix3D.size(); i++){
            multiplyMatrix(multiplicator, matrix3D.get(i));
        }
    }

    public static void multiplyMatrix(double multiplicator, List<List<Double>> matrix) {
        for(int i=0; i<matrix.size(); i++){
            multiplyVector(multiplicator, matrix.get(i));
        }
    }

    public static void multiplyVector(double multiplicator, List<Double> vector) {
        for(int i=0; i<vector.size(); i++){
            vector.set(i, vector.get(i)*multiplicator);
        }
    }

    public static List<List<List<Double>>> create3DMatrix(int sizeX, int sizeY, int sizeZ, double filler) {
        List<List<List<Double>>> matrix3D = new ArrayList();
        for(int i=0; i<sizeZ; i++){
            matrix3D.add(createMatrix(sizeX,sizeY,filler));
        }
        return matrix3D;
    }

    public static List<List<Double>> createMatrix(int sizeX, int sizeY, double filler){
        List<List<Double>> matrix = new ArrayList();
        for(int line=0; line<sizeY; line++){ // Add sizeY rows (to have sizeY columns)
            matrix.add(createVector(sizeX,filler));
        }
        return matrix;
    }

    public static List<Double> createVector(int size, double filler) {
        List<Double> vector = new ArrayList();
        for(int line=0; line<size; line++){
            vector.add(filler);
        }
        return vector;
    }

    public static List<Double> createVectorLoop(double startValue, double increaseAmount, double endValue){
        List<Double> resultArray = new ArrayList();
        double currentValue = startValue;
        while(currentValue <= endValue){
            resultArray.add(currentValue);
            currentValue+=increaseAmount;
        }
        return resultArray;
    }

    public static double maxVector(List<Double> vector) {
        double max = vector.get(0);
        for (Double d :
                vector) {
            max = Math.max(d,max);
        }
        return max;
    }

    public static double minVector(List<Double> vector) {
        double min = vector.get(0);
        for (Double d :
                vector) {
            min = Math.min(d,min);
        }
        return min;
    }

    public static void setMatrix3DValue(List<List<List<Double>>> matrix3D, int row, int col, int page, double value) {
        matrix3D.get(page).get(col).set(row,value);
    }

    public static List<List<Double>> getPage(List<List<List<Double>>> matrix3D, int pageIndex) {
        return  matrix3D.get(pageIndex);
    }

    /**
     * Create an under matrix of columns from startIndex to stopIndex containing all rows
     * @param matrix
     * @param startIndex
     * @param stopIndex
     * @return
     */
    public static List<List<Double>> getUnderMatrixColumnIteration(List<List<Double>> matrix, int startIndex, int stopIndex) {
        List<List<Double>> resultMatrix = new ArrayList<>();
        for(int col=startIndex; col <=stopIndex; col++){
            resultMatrix.add(getColumn(matrix, col));
        }
        return resultMatrix;
    }

    public static List<Double> getColumn(List<List<Double>> matrix, int index) {
        List<Double> col = createVector(matrix.size(), 0);
        for(int r=0; r < matrix.size(); r++){
            col.set(r,matrix.get(r).get(index));
        }
        return col;
    }

    private static List<Double> getRow(List<List<Double>> matrix, int index) {
        return matrix.get(index);
    }

    public static double matrixSum(List<List<Double>> matrix) {
        double sum = 0;
        for (List<Double> row :
                matrix) {
            sum += vectorSum(row);
        }
        return sum;
    }

    private static double vectorSum(List<Double> vector) {
        double sum = 0;
        for (Double value :
                vector) {
            sum += value;
        }
        return sum;
    }

    public static List<Double> cumsum(List<Double> vector) {
        double sum = 0;
        List<Double> resultVector = vector;
        for (int i=0; i<resultVector.size(); i++) {
            sum+= resultVector.get(i);
            resultVector.set(i,sum);
        }
        return resultVector;
    }

    /**
     * Create a vector of 0 and 1.
     * 0 at i if the element at position i is smaller than valueToTest
     * @param vector
     * @param valueToTest
     * @return
     */
    public static List<Double> compareBiggerEqual(List<Double> vector, double valueToTest) {
        List<Double> resultVector = vector;
        for (int i=0; i<resultVector.size(); i++) {
            resultVector.set(i,(resultVector.get(i)>=valueToTest? 1. : 0.));
        }
        return resultVector;
    }

    public static List<Double> getIndexOfNonZeros(List<Double> vector) {
        List<Double> resultVector = new ArrayList<>();
        for (int i=0; i<vector.size(); i++) {
            if(vector.get(i)!=0)
                resultVector.add((double) i);
        }
        return resultVector;
    }
}
