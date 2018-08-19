import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVPrinter;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;

public class CSVWriter {
    private String outputPath = "./sample.csv";

    public void CSVWriter(String path){
        outputPath = path;
    }

    public void CSVWriter(){
    }

    public void write() throws IOException {
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
        }
    }

}
