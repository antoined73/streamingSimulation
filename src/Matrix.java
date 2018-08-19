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

    public static double vectorSum(List<Double> vector) {
        double sum = 0;
        for (Double value :
                vector) {
            sum += value;
        }
        return sum;
    }

    public static List<Double> cumsum(List<Double> vector) {
        double sum = 0;
        List<Double> resultVector = new ArrayList<>(vector);
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
        List<Double> resultVector = new ArrayList<>(vector);
        for (int i=0; i<resultVector.size(); i++) {
            resultVector.set(i,(resultVector.get(i)>=valueToTest? 1. : 0.));
        }
        return resultVector;
    }

    /**
     * Create a vector of 0 and 1.
     * 0 at i if the element at position i is equal than valueToTest
     * @param vector
     * @param valueToTest
     * @return
     */
    public static List<Double> compareEqual(List<Double> vector, double valueToTest) {
        List<Double> resultVector = new ArrayList<>(vector);
        for (int i=0; i<resultVector.size(); i++) {
            resultVector.set(i,(resultVector.get(i)>=valueToTest? 1. : 0.));
        }
        return resultVector;
    }

    public static List<Integer> getIndexOfNonZeros(List<Double> vector) {
        List<Integer> resultVector = new ArrayList<>();
        for (int i=0; i<vector.size(); i++) {
            if(vector.get(i)!=0)
                resultVector.add(i);
        }
        return resultVector;
    }

    public static List<Double> randomVector(int size) {
        List<Double> vector = createVector(size,0);
        for (int i=0; i<vector.size(); i++){
            vector.set(i, Math.random());
        }
        return vector;
    }

    public static void addValues(List<Double> vectorResult, List<Double> addedVector) {
        for (int i=0; i<vectorResult.size(); i++){
            vectorResult.set(i, vectorResult.get(i) + addedVector.get(i));
        }
    }

    public static List<List<Double>>  repeatMatrix(List<Double> matrix, int nbCopiesOnY, int nbCopiesOnX) {
        List<List<Double>> resultMatrix = createMatrix(0,0,0);
        resultMatrix.add(cloneVector(matrix));

        for (int x=0; x<nbCopiesOnX; x++){
            for (int i=0; i<resultMatrix.size(); i++){
                resultMatrix.get(i).addAll(cloneVector(resultMatrix.get(i)));
            }
        }

        List<List<Double>> tempMatrix = cloneMatrix(resultMatrix);
        for (int y=0; y<nbCopiesOnY; y++){
            resultMatrix.addAll(cloneMatrix(tempMatrix));
        }
        return resultMatrix;
    }

    private static List<List<Double>> cloneMatrix(List<List<Double>> matrix) {
        List<List<Double>> clonedMatrix = new ArrayList<>();
        for (List<Double> row:
             matrix) {
            clonedMatrix.add(cloneVector(row));
        }
        return clonedMatrix;
    }

    public static List<Double> cloneVector(List<Double> vector) {
        List<Double> clonedVector = new ArrayList<>();
        for (Double d :
                vector) {
            clonedVector.add(d);
        }
        return clonedVector;
    }

    public static List<Double> cloneVectorToDouble(List<Integer> vector) {
        List<Double> clonedVector = new ArrayList<>();
        for (Integer d :
                vector) {
            clonedVector.add((double)d);
        }
        return clonedVector;
    }

    public static void substractValuesMatrix(List<List<Double> > resultMatrix, List<List<Double>> substractedMatrix) {
        for (int i=0; i<resultMatrix.size(); i++){
            substractValuesVector(resultMatrix.get(i), substractedMatrix.get(i));
        }
    }

    private static void substractValuesVector(List<Double> resultVector, List<Double> substractedVector) {
        for (int i=0; i<resultVector.size(); i++){
            resultVector.set(i, resultVector.get(i) - substractedVector.get(i));
        }
    }

    public static void applyPowerMatrix(List<List<Double>> matrix, int power) {
        for (int i=0; i<matrix.size(); i++){
            applyPowerVector(matrix.get(i), power);
        }
    }

    private static void applyPowerVector(List<Double> vector, int power) {
        for (int i=0; i<vector.size(); i++){
            vector.set(i, Math.pow(vector.get(i),power));
        }
    }

    public static List<Double> sumOfEachRow(List<List<Double>> matrix) {
        List<Double> vector = createVector(matrix.size(),0);
        for (int i=0; i<vector.size(); i++){
            vector.set(i, vectorSum(matrix.get(i)));
        }
        return vector;
    }

    public static List<Double> substract(double value, List<Double> vector) {
        List<Double> resultVector = cloneVector(vector);
        for (int i=0; i<vector.size(); i++){
            resultVector.set(i, value-resultVector.get(i));
        }
        return resultVector;
    }

    /**
     * add value to each item of vector
     * @param value
     * @param vector
     * @return
     */
    public static List<Double> add(double value, List<Double> vector) {
        List<Double> resultVector = cloneVector(vector);
        for (int i=0; i<vector.size(); i++){
            resultVector.set(i, value+resultVector.get(i));
        }
        return resultVector;
    }

    public static List<Double> divideMatrix(List<Double> vector, double divider) {
        List<Double> resultVector = cloneVector(vector);
        for(int i=0; i<resultVector.size(); i++){
            resultVector.set(i, resultVector.get(i)/divider);
        }
        return resultVector;
    }

    public static List<List<Double>> minimalMatrix(List<List<List<Double>>> matrix3D) {
        List<List<Double>>  resultMatrix = cloneMatrix(matrix3D.get(0));
        for (List<List<Double>> matrix :
                matrix3D) {
            if(compareBiggerEqual(resultMatrix, matrix)){
                resultMatrix = matrix;
            }
        }
        return resultMatrix;
    }

    /**
     * Return true if resultMatrix is Bigger or equal than matrix
     * @param resultMatrix
     * @param matrix
     * @return
     */
    private static boolean compareBiggerEqual(List<List<Double>> resultMatrix, List<List<Double>> matrix) {
        return matrixSum(resultMatrix) >= matrixSum(matrix);
    }

    private static Double minMatrix(List<List<Double>> matrix) {
        double min = matrix.get(0).get(0);
        for (List<Double> vector :
                matrix) {
            Math.min(min, minVector(vector));
        }
        return min;
    }

    public static List<Double> minValueOfEachRow(List<List<Double>> matrix) {
        List<Double> resultVector = new ArrayList<>();
        for (List<Double> v :
                matrix) {
            resultVector.add(minVector(v));
        }
        return resultVector;
    }

    /**
     * Keep only the values of the vector that are in front of 1 in the filter.
     * @param vector
     * @param filter
     * @return
     */
    public static List<Integer> applyBooleanFilter(List<Double> vector, List<Double> filter) {
        List<Integer> resultVector = new ArrayList<>();
        for (int i=0; i<vector.size(); i++) {
            if(filter.get(i) != 0.){
                resultVector.add(i);
            }
        }
        return resultVector;
    }

    public static void setValueAtIndexes(List<Double> vector, List<Integer> indexes, double value) {
        for (int i=0; i<indexes.size(); i++){
            vector.set(indexes.get(i),value);
        }
    }

    public static List<List<Double>> multiplyElementByElementMatrix(List<List<Double>> m1, List<List<Double>> m2) {
        List<List<Double>> resultMatrix = cloneMatrix(m1);
        for (int i=0; i<resultMatrix.size(); i++){
            if(i<m1.size()&& i< m2.size())
                resultMatrix.set(i, multiplyElementByElementVector(m1.get(i),m2.get(i)) );
        }
        return resultMatrix;
    }

    private static List<Double> multiplyElementByElementVector(List<Double> v1, List<Double> v2) {
        List<Double> resultVector = cloneVector(v1);
        for (int i=0; i<resultVector.size(); i++){
            resultVector.set(i, v1.get(i) * v2.get(i) );
        }
        return resultVector;
    }

    public static void setValueInMatrix3DToElementsSmallerThan(List<List<List<Double>>> matrix3D, double compareValue, double value) {
        for (List<List<Double>> matrix :
                matrix3D) {
            setValueInMatrixToElementsSmallerThan(matrix, compareValue, value);
        }
    }

    private static void setValueInMatrixToElementsSmallerThan(List<List<Double>> matrix, double compareValue, double value) {
        for (List<Double> vector :
                matrix) {
            setValueInVectorToElementsSmallerThan(vector, compareValue, value);
        }
    }

    private static void setValueInVectorToElementsSmallerThan(List<Double> vector, double compareValue, double value) {
        for (int i=0; i<vector.size(); i++) {
            if(vector.get(i) < compareValue){
                vector.set(i,value);
            }
        }
    }

    public static void setValueInMatrix3DToElementsEqualTo(List<List<List<Double>>> matrix3D, double compareValue, double value) {
        for (List<List<Double>> matrix :
                matrix3D) {
            setValueInMatrixToElementsEqualTo(matrix, compareValue, value);
        }
    }

    private static void setValueInMatrixToElementsEqualTo(List<List<Double>> matrix, double compareValue, double value) {
        for (List<Double> vector :
                matrix) {
            setValueInVectorToElementsEqualTo(vector, compareValue, value);
        }
    }

    private static void setValueInVectorToElementsEqualTo(List<Double> vector, double compareValue, double value) {
        for (int i=0; i<vector.size(); i++) {
            if(vector.get(i) == compareValue){
                vector.set(i,value);
            }
        }
    }

    public static List<List<Double>> sumOfEachMatrix(List<List<List<Double>>> matrix3D) {
        List<List<Double>> resultMatrix = createMatrix(matrix3D.get(0).get(0).size(),matrix3D.get(0).size(),0);
        for (List<List<Double>> matrix :
                matrix3D) {
            for (int y=0; y<matrix.size(); y++) {
                for (int x=0; x<matrix.get(y).size(); x++) {
                    resultMatrix.get(y).set(x,
                            resultMatrix.get(y).get(x) + matrix.get(y).get(x));
                }
            }
        }
        return  resultMatrix;
    }

    public static List<List<List<Double>>> sum3DOfEachRowInMatrix3D(List<List<List<Double>>> matrix3D) {
        List<List<List<Double>>> resultMatrix3D = new ArrayList<>();
        for (List<List<Double>> matrix :
             matrix3D) {
            List<List<Double>> m = new ArrayList<>();
            m.add(sumOfEachRow(matrix));
            resultMatrix3D.add(m);
        }
        return resultMatrix3D;
    }

    public static List<Double> vectorSumOfEachRowInMatrix(List<List<Double>> matrix) {
        List<Double> resultVector = new ArrayList<>();
        for (int line=0; line<matrix.size(); line++) {
            resultVector.add(Matrix.vectorSum(matrix.get(line)));
        }
        return resultVector;
    }

    public static List<Double> getValuesSmallerThan(List<Double> vector, double compareValue) {
        List<Double> resultVector = new ArrayList<>();
        for (int i=0; i<vector.size(); i++) {
            if(vector.get(i) < compareValue){
                resultVector.add(vector.get(i));
            }
        }
        return resultVector;
    }

    public static double matrix3Dsum(List<List<List<Double>>> matrix3D) {
        double sum = 0;
        for (List<List<Double>> matrix : matrix3D){
            sum+= matrixSum(matrix);
        }
        return sum;
    }

    public static List<List<List<Double>>> multiplyElementByElementMatrix3D(List<List<List<Double>>> matrix3D_1, List<List<List<Double>>> matrix3D_2) {
        List<List<List<Double>>> resultMatrix3D = new ArrayList<>();
        for (int i=0; i<matrix3D_1.size(); i++) {
            if(i < matrix3D_1.size() && i < matrix3D_2.size())
                resultMatrix3D.add(multiplyElementByElementMatrix(matrix3D_1.get(i),matrix3D_2.get(i)));
        }
        return resultMatrix3D;
    }

    public static List<List<List<Integer>>> create3DMatrixInteger(int sizeX, int sizeY, int sizeZ, int filler) {
        List<List<List<Integer>>> matrix3D = new ArrayList();
        for(int i=0; i<sizeZ; i++){
            matrix3D.add(createMatrixInteger(sizeX,sizeY,filler));
        }
        return matrix3D;
    }

    public static List<List<Integer>> createMatrixInteger(int sizeX, int sizeY, int filler){
        List<List<Integer>> matrix = new ArrayList();
        for(int line=0; line<sizeY; line++){ // Add sizeY rows (to have sizeY columns)
            matrix.add(createVectorInteger(sizeX,filler));
        }
        return matrix;
    }

    public static List<Integer> createVectorInteger(int size, int filler) {
        List<Integer> vector = new ArrayList();
        for(int line=0; line<size; line++){
            vector.add(filler);
        }
        return vector;
    }

    public static List<List<List<Integer>>> DoubleMatrix3DtoInteger(List<List<List<Double>>> m) {
        List<List<List<Integer>>> matrix3D = new ArrayList();
        for(int i=0; i<m.size(); i++){
            matrix3D.add(cloneMatrixInteger(m.get(i)));
        }
        return matrix3D;
    }

    public static List<List<Integer>> cloneMatrixInteger(List<List<Double>> m) {
        List<List<Integer>> matrix = new ArrayList();
        for(int line=0; line<m.size(); line++){ // Add sizeY rows (to have sizeY columns)
            matrix.add(cloneVectorToInteger(m.get(line)));
        }
        return matrix;
    }

    public static List<Integer> cloneVectorToInteger(List<Double> v) {
        List<Integer> vector = new ArrayList();
        for(int line=0; line<v.size(); line++){
            vector.add(v.get(line).intValue());
        }
        return vector;
    }

    public static List<Integer> cloneVectorInteger(List<Integer> v) {
        List<Integer> vector = new ArrayList();
        for(int line=0; line<v.size(); line++){
            vector.add(v.get(line));
        }
        return vector;
    }

    public static void printMatrix3D(List<List<List<Double>>> matrix3D) {
        for(int i=0; i<matrix3D.size(); i++){
            System.out.println("feuille "+i);
            printMatrix(matrix3D.get(i));
        }
    }

    private static void printMatrix(List<List<Double>> matrix) {
        for(int i=0; i<matrix.size(); i++){
            printVector(matrix.get(i));
        }
    }

    public static void printVector(List<Double> vector) {
        String line="[ ";
        for(int i=0; i<vector.size(); i++){
            line+= vector.get(i);
            if(i+1 < vector.size())
                line+= ", ";
        }
        line+= " ]";
        System.out.println(line);
    }

    public static List<Integer> getIndexOfValuesSmallerThan(List<Double> vector, double value) {
        List<Integer> resultIndexes = new ArrayList<>();
        for (int i=0; i<vector.size(); i++) {
            if(vector.get(i) < value)
                resultIndexes.add(i);
        }
        return resultIndexes;
    }

    public static List<List<Double>> maximalMatrix(List<List<List<Double>>> matrix3D) {
        List<List<Double>>  resultMatrix = cloneMatrix(matrix3D.get(0));
        for (List<List<Double>> matrix :
                matrix3D) {
            if(!compareBiggerEqual(resultMatrix, matrix)){
                resultMatrix = matrix;
            }
        }
        return resultMatrix;
    }

    public static void setValueInMatrix3DToElementsBiggerThan(List<List<List<Double>>> matrix3D, int compareValue, double value) {
        for (List<List<Double>> matrix :
                matrix3D) {
            setValueInMatrixToElementsBiggerThan(matrix, compareValue, value);
        }
    }

    private static void setValueInMatrixToElementsBiggerThan(List<List<Double>> matrix, int compareValue, double value) {
        for (List<Double> vector :
                matrix) {
            setValueInVectorToElementsBiggerThan(vector, compareValue, value);
        }
    }

    private static void setValueInVectorToElementsBiggerThan(List<Double> vector, int compareValue, double value) {
        for (int i=0; i<vector.size(); i++) {
            if(vector.get(i)>compareValue){
                vector.set(i, value);
            }
        }
    }
}
