import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVPrinter;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class CSVWriter {
    private String outputPath = "./sample.csv";
    private String inputPath = "D:/invasion/";
    public int maxSegmentNumber = 0;
    public List<String> qualities;
    public int maxX = 0;
    public int maxY = 0;

    public CSVWriter(String outputPath, String m4sDirectory){
        this.outputPath = outputPath;
        this.inputPath = m4sDirectory;
        qualities = new ArrayList<>();
    }

    public CSVWriter(){
        qualities = new ArrayList<>();
    }


    public void write() throws IOException {

        try (
                BufferedWriter writer = Files.newBufferedWriter(Paths.get(outputPath));

                CSVPrinter csvPrinter = new CSVPrinter(writer, CSVFormat.DEFAULT
                        .withHeader("Quality", "Segment", "X", "Y", "Size"));
        ) {
            File folder = new File(inputPath);
            File[] listOfFiles = folder.listFiles();

            for (int i = 0; i < listOfFiles.length; i++) {
                if (listOfFiles[i].isFile() && listOfFiles[i].getName().endsWith(".m4s") && !listOfFiles[i].getName().contains("audio")) {
                    //invasion-720x576-0111__track1_168.m4s

                    String fileName = listOfFiles[i].getName().substring(0,listOfFiles[i].getName().length()-4);
                    //quality
                    String[] fileNameParts = fileName.split("-");
                    if(!qualities.contains(fileNameParts[1])){
                        System.out.println(fileNameParts[1]);
                        if(qualities.size()==0)
                            qualities.add(fileNameParts[1]);
                        else{
                            String[] res = fileNameParts[1].split("x");
                            int current_w = Integer.parseInt(res[0]);
                            int current_h = Integer.parseInt(res[1]);

                            boolean added = false;
                            for (int qi=0; qi<qualities.size(); qi++) {
                                String[] resolutions = qualities.get(qi).split("x");
                                int w = Integer.parseInt(resolutions[0]);
                                int h = Integer.parseInt(resolutions[1]);
                                //System.out.println(w+","+h);
                                if(current_h*current_w <= h*w){
                                    qualities.add(qi,fileNameParts[1]);
                                    added = true;
                                    break;
                                }
                            }

                            if(!added){
                                qualities.add(fileNameParts[1]);
                            }
                        }
                    }
                }
            }

            for (int i = 0; i < listOfFiles.length; i++) {
                if (listOfFiles[i].isFile() && listOfFiles[i].getName().endsWith(".m4s") && !listOfFiles[i].getName().contains("audio")) {
                    //invasion-720x576-0111__track1_168.m4s

                    String fileName = listOfFiles[i].getName().substring(0,listOfFiles[i].getName().length()-4);

                    //quality
                    String[] fileNameParts = fileName.split("-");
                    int q = 0;

                    for(int qi=0; qi<qualities.size(); qi++){ //get current quality index
                        if(qualities.get(qi).equals(fileNameParts[1])){
                            q = qi;
                        }
                    }

                    //x, y & segment
                    int segNb, x, y;
                    String[] lastParameters = fileNameParts[2].split("_");
                    segNb = Integer.parseInt(lastParameters[3]);
                    x = Integer.parseInt(String.valueOf(lastParameters[0].charAt(0)));
                    y = Integer.parseInt(String.valueOf(lastParameters[0].charAt(1)));

                    maxSegmentNumber = Math.max(segNb, maxSegmentNumber);
                    maxX = Math.max(x, maxX);
                    maxY = Math.max(y, maxY);
                    //size
                    int size = (int) listOfFiles[i].length();

                    csvPrinter.printRecord(q, segNb, x, y, size);
                }
            }

            csvPrinter.flush();
        }



        /***
        try (
                BufferedWriter writer = Files.newBufferedWriter(Paths.get(outputPath));

                CSVPrinter csvPrinter = new CSVPrinter(writer, CSVFormat.DEFAULT
                        .withHeader("Quality", "Segment", "X", "Y", "Size"));
        ) {
            for (int q=0; q<2; q++){
                for (int s=0; s<20; s++){
                    for (int x=0; x<2; x++){
                        for (int y=0; y<2; y++){
                            csvPrinter.printRecord(q, s, x, y,(1+q)*150000);
                        }
                    }
                }
            }

            csvPrinter.flush();
        }***/
    }

}
