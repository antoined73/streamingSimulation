import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;

import java.io.IOException;
import java.io.Reader;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

public class CSVReader {
    private String filePath = "./sample.csv";

    private int nb_tiles=4;
    private int nb_segments=20;
    private int nb_levels=2;

    public CSVReader(String filePath, int nbtiles, int nbsegs, int nblvls){
        this.nb_segments = nbsegs;
        this.nb_tiles = nbtiles;
        this.nb_levels = nblvls;
        this.filePath = filePath;
    }

    public CSVReader(String filePath){
        this.filePath = filePath;
    }

    public List<List<List<Double>>> read() throws IOException {
        List<List<List<Double>>> matrix3D = Matrix.create3DMatrix(nb_segments,nb_tiles,nb_levels,0);
        try (
                Reader reader = Files.newBufferedReader(Paths.get(filePath));
                CSVParser csvParser = new CSVParser(reader, CSVFormat.DEFAULT
                        .withFirstRecordAsHeader()
                        .withIgnoreHeaderCase()
                        .withTrim());
        ) {
            for (CSVRecord csvRecord : csvParser) {
                // Accessing values by Header names
                int quality = Integer.parseInt(csvRecord.get("Quality"));
                int segment = Integer.parseInt(csvRecord.get("Segment"));
                int x = Integer.parseInt(csvRecord.get("X"));
                int y = Integer.parseInt(csvRecord.get("Y"));
                int size = Integer.parseInt(csvRecord.get("Size"));


                System.out.println("Quality : " + quality+" Segment : " + segment+" X : " + x+" Y : " + y+ " Size : " + size);

                fillMatrix3D(matrix3D, quality, segment-1, x, y, size);
            }
        }
        return matrix3D;
    }

    List<String> indexs_xy =  new ArrayList<>();
    private void fillMatrix3D(List<List<List<Double>>> matrix3D, int quality, int segment, int x, int y, int size) {

        List<List<Double>> m = null;
        if(quality >= matrix3D.size()){
            System.err.println("quality greater than given in initialisation : "+quality);
        }else {
            m = matrix3D.get(quality);
        }

        int indexTile = getTileIndex(indexs_xy, m, x, y);
        //System.out.println(x+","+y+": "+indexTile);
        List<Double> line = null;
        if(indexTile >= m.size()){ //add line for tile
            System.err.println("number of tiles greater than given in initialisation : "+x+","+y);
        }else{
            line = m.get(indexTile);
        }

        if(segment>line.size()){
            System.err.println("number of segments greater than given in initialisation : "+nb_segments);
        }else{
            line.set(segment, (double) size);
        }
    }

    private int getTileIndex(List<String> indexs_xy, List<List<Double>> m, int x, int y) {
        for(int i=0; i<indexs_xy.size(); i++){
            if(indexs_xy.get(i).equals(x + "" + y))
                return i;
        }
        indexs_xy.add(x+""+y);
        return indexs_xy.size()-1;
    }


}