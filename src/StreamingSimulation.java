import java.util.*;

import static java.lang.Math.floor;
import static java.lang.Math.sqrt;

public class StreamingSimulation {

    int nb_of_levels = 2;
    private int K_lookahead = 3;
    int nb_of_tiles = 4;
    int nb_of_segments = 20;
    double seg_duration = 1;
    double deltaDownload = seg_duration * K_lookahead;

    List<List<Double>> tile_centers_xy; // matrix of 4 lines containing each 2 columns
    double tilecenter_offset;
    List<Double> x, y; //float[2] sizes of the video square


    public StreamingSimulation(){
        tile_centers_xy = Matrix.createMatrix(2,nb_of_tiles,0);
        tilecenter_offset = 1/(2*sqrt(nb_of_tiles));

        // assuming sphere as 1X1 square
        x = Matrix.createVectorLoop(tilecenter_offset, 1/sqrt(nb_of_tiles), 1-tilecenter_offset);
        y = x;

        int currentTileNb = 1;
        for (int ix=0; ix < x.size(); ix++){
            for (int iy=0; iy < y.size(); iy++){
                List<Double> line = new ArrayList();
                line.add(x.get(ix));
                line.add(y.get(iy));

                tile_centers_xy.set(currentTileNb-1,line);
                currentTileNb++;
            }
        }
    }

    public void startStreaming(){

        double tile_size = 300 * Math.pow(10,3); //TODO : Get Tile size in .csv
        double duration= 10 * nb_of_segments * seg_duration;
        double cst_bw = 1* Math.pow(10,6);
        int Bmin=2;
        int Bmax=5;

        List<List<List<Double>>> s_ijl = Matrix.create3DMatrix(nb_of_segments,nb_of_tiles,nb_of_levels,1);
        Matrix.multiply3D(tile_size, s_ijl);

        List<List<Double>> tmp_matrix = Matrix.getPage(s_ijl,2-1);
        Matrix.multiplyMatrix(0.5, tmp_matrix );
        s_ijl.set(1,tmp_matrix);


        List<Double> bw_trace = Matrix.createVector(((int) floor(duration/seg_duration)), 1);
        Matrix.multiplyVector(cst_bw,bw_trace);

        //---Streaming start---//
        double time_user=0; // counted in number of segments
        double time_video=0; // counted in number of segments played out so far
        int seg_played=0;

        List<List<List<Double>>> buf_it = Matrix.create3DMatrix(Bmax+2,nb_of_tiles,nb_of_levels,1);
        Matrix.multiply3D(-3, buf_it); // stores the segments index, j, 1st level 0 means no seg
        int nb_of_startupseg = (int) floor(Bmin/2f);

        List<List<Double>> page = Matrix.getPage(buf_it, 0);
        for(int col=0; col < nb_of_startupseg; col++){ //iterate throught the columns from 0 to nb_of_startupseg
            int currentNumber = 1;
            List<Double> column = Matrix.getColumn(page,col);
            for(int row=0; row < column.size(); row++){
                Matrix.setMatrix3DValue(buf_it, col, row,0, currentNumber);
                //currentNumber++;
            }
        }
        // Startup period
        List<List<Double>> played_qualities = Matrix.createMatrix(nb_of_tiles,nb_of_segments,0);
        List<List<Double>> dlded_size = Matrix.getUnderMatrixColumnIteration( Matrix.getPage(s_ijl,1), 1-1 , nb_of_startupseg-1);
        double dlded_size_sum = Matrix.matrixSum(dlded_size);

        /**Works from beginning of file to here**/

        time_user = computeDlDelay(dlded_size_sum,bw_trace,0); // startup delay

        //-- update bw estimate
        List<Double> truedl_rate = bw_trace.subList(0, (int) time_user);
        double Ct = estimateNextBtw(1,truedl_rate);

        List<Double> FoV_xy = Matrix.createVector(2,0);
        double total_stalltime=0;


        while (seg_played<nb_of_segments){
            //-- prepare for dl decision
            FoV_xy = generateNewFOV(FoV_xy);
            List<List<Double>> repmat = Matrix.repeatMatrix(FoV_xy,nb_of_tiles-1,0);
            Matrix.substractValuesMatrix(repmat,tile_centers_xy);
            Matrix.applyPowerMatrix(repmat,2);
            List<Double> dist_FoVtiles = Matrix.sumOfEachRow(repmat);

            List<Double> diff = Matrix.substract(Matrix.maxVector(dist_FoVtiles),dist_FoVtiles);
            List<Double> p_ij = Matrix.divideMatrix(diff,Matrix.vectorSum(diff));


            List<List<Double>> j_ti = Matrix.minimalMatrix(buf_it);
            List<Double> j_ti_min = Matrix.minValueOfEachRow(j_ti);

            //j_ti(j_ti==(nb_of_segments+3))=time_video;
            List<Double> boolVector = Matrix.compareEqual(j_ti_min,(nb_of_segments+3));
            List<Double> indexes = Matrix.applyBooleanFilter(j_ti_min, boolVector);
            Matrix.setValueAtIndexes(j_ti_min,indexes,time_video);

            Matrix.add(1,j_ti_min);
            double j_t = Matrix.minVector(j_ti_min);


            //-- make dl decision
            List<List<List<Double>>> x_ijl = instant_optim(K_lookahead,deltaDownload,p_ij,s_ijl,Ct,buf_it,j_ti_min,Bmin,Bmax);

            //-- buffers states after dl finished (right before next download attempt)

            dlded_size = Matrix.multiplyElementByElementMatrix(x_ijl.get(0),s_ijl.get(0));

            dlded_size_sum = Matrix.matrixSum(dlded_size);
            double time_to_dl = computeDlDelay(dlded_size_sum,bw_trace,time_user); //counted in segments

            //-- temporary transform
            List<List<List<Double>>> buf_tmp = buf_it;
            Matrix.setValueInMatrix3DToElementsSmallerThan(buf_tmp,(nb_of_segments+3),1);
            Matrix.setValueInMatrix3DToElementsEqualTo(buf_tmp,(nb_of_segments+3),0);

            //-- end of temporary transform

            double cursizebuf_min= Matrix.minVector(Matrix.sumOfEachRow(Matrix.sumOfEachMatrix(buf_tmp)));
            double stall_time = Math.max(0,time_to_dl-cursizebuf_min);
            total_stalltime=total_stalltime+stall_time;


            List<List<List<Integer>>> buf_it_int = Matrix.DoubleMatrix3DtoInteger(buf_it);
            List<List<Integer>> played_qualities_int = Matrix.cloneMatrixInteger(played_qualities);
            List<List<List<Integer>>> x_ijl_int =  Matrix.DoubleMatrix3DtoInteger(x_ijl);
            List<Integer> j_ti_min_int = Matrix.cloneVectorInteger(j_ti_min);
            updateBufferStates(buf_it_int,played_qualities_int,x_ijl_int,j_ti_min_int,time_to_dl,K_lookahead,time_video,Bmax);

            //-- update bw estimate
            truedl_rate = bw_trace.subList((int)time_user-1, ((int)(time_user+time_to_dl-1-1)) );
            Ct = estimateNextBtw(Ct,truedl_rate);

            //-- update times
            time_user = time_user + deltaDownload + Math.ceil(stall_time/seg_duration);
            time_video = time_video+time_to_dl-stall_time;

        }

    }


    private void updateBufferStates(List<List<List<Integer>>> buf_it, List<List<Integer>> played_qualities, List<List<List<Integer>>> x_ijl, List<Integer> j_ti,
                                    double time_to_dl, int K_lookahead, double time_video, int Bmax){ //update buf_it and played_qualities
        this.nb_of_tiles = x_ijl.size();
        this.nb_of_segments = x_ijl.get(0).size();

        int j_t = j_ti.stream().min(Comparator.naturalOrder()).orElse(0);

        //temporary transform
        List<List<List<Integer>>> buf_tmp = new ArrayList<>(buf_it);
        for (List<List<Integer>> aBuf_tmp1 : buf_tmp) {
            for (int j = 0; j < aBuf_tmp1.size(); ++j) {
                for (int k = 0; k < aBuf_tmp1.get(j).size(); ++k) {
                    if (aBuf_tmp1.get(j).get(k) < nb_of_segments + 3) {
                        aBuf_tmp1.get(j).set(k, 1);
                    } else if (aBuf_tmp1.get(j).get(k) == nb_of_segments + 3) {
                        aBuf_tmp1.get(j).set(k, 0);
                    }
                }
            }
        }
        //end of temporary transform

        List<Integer> bufsize_i = new ArrayList<>();
        {
            List<List<Integer>> tmp = new LinkedList<>();
            for (List<List<Integer>> aBuf_tmp : buf_tmp) {
                List<Integer> tmp2 = new LinkedList<>();
                for (List<Integer> anABuf_tmp : aBuf_tmp) {
                    tmp2.add(anABuf_tmp.stream().mapToInt(t -> t).sum());
                }
                tmp.add(tmp2);
            }
            for (List<Integer> aTmp : tmp) {
                bufsize_i.add(aTmp.stream().mapToInt(t -> t).sum());
            }
        }

        //-- first fill: everything that was scheduled for dl
        buf_tmp = new ArrayList<>(buf_it);
        for (int i = 0 ; i < nb_of_tiles; ++i) {
            int ind_placeinbuf = 0;
            for (int j = j_ti.get(i); j < j_t + K_lookahead -1; ++j) {
            	for (int k = 0 ;
                     k < x_ijl.get(i).get(j).size() ; ++k) {
            	    if (x_ijl.get(i).get(j).get(k) == 1) {
            	        buf_tmp.get(i).get(bufsize_i.get(i) + ind_placeinbuf).set(k, j);
                    }
                }
                ind_placeinbuf++;
            }
        }

        //-- then drain: everything that has been read in max(deltaDownload,time_to_dl)
	    for (int played_seg = (int)time_video ; played_seg < (int)time_video + time_to_dl ; ++played_seg) {
            for (int i = 0 ; i < nb_of_tiles ; ++i) {
                int finalPlayed_seg = played_seg;
                int finalI = i;
                buf_tmp.get(i).get(played_seg - (int)time_video)
                        .stream()
                        .filter(d -> d.equals(finalPlayed_seg))
                        .findFirst()
                        .ifPresent(t -> played_qualities.get(finalI).set(finalPlayed_seg, t));
            }
        }

        buf_it.clear();
	    for (int i = 0 ; i < buf_tmp.size() ; ++i) {
            List<List<Integer>> tmp = new ArrayList<>();
	        int lim = (int) Math.min(time_to_dl + Bmax + 2, buf_tmp.get(i).size());
	    	for (int j = (int) time_to_dl; j < lim ; ++j) {
				tmp.add(new ArrayList<>(buf_tmp.get(i).get(j)));
            }
            buf_it.add(tmp);
        }
        for (int i = 0 ; i < buf_it.size() ; ++i) {
            for (int j = buf_it.get(i).size() ; j < Bmax + 2 ; ++j) {
                for (int k = 0 ; k < buf_it.get(i).get(j).size() ; ++k) {
                    buf_it.get(i).get(j).set(k, nb_of_segments + 3);
                }
            }
        }
    }

    private List<Double> generateNewFOV(List<Double> FoV_xy){
        List<Double> newFov = new ArrayList<>(FoV_xy);
        Matrix.addValues(newFov, Matrix.randomVector(2));
        //FoV_xy(FoV_xy>1) = FoV_xy(FoV_xy>1)-1;
        for (int i=0; i<newFov.size(); i++){
            if(newFov.get(i)>1)
                newFov.set(i,newFov.get(i)-1);
        }
        return newFov;
    }

    private double estimateNextBtw(double Ct, List<Double> truedl_rate) { //update Ct
        double newCt = Ct;

        float alpha = 1;
        for(int t=0; t<truedl_rate.size(); t++){
                newCt = ((1-alpha)* newCt) + (alpha * truedl_rate.get(t));
        }
        return newCt;
    }

    private double computeDlDelay(double dlded_size,List<Double> bw_trace,double start_time_user){ //update time_user
        List<Double> trace = new ArrayList<>(bw_trace.subList((int) start_time_user,bw_trace.size()));

        List<Double> size_dldable = Matrix.cumsum(trace);
        List<Double> tmp = Matrix.compareBiggerEqual(size_dldable,dlded_size);
        List<Double> time_dlded = Matrix.getIndexOfNonZeros(tmp);

        double time = (time_dlded.size() > 0 ? time_dlded.get(0) : 0);

        return time+1;
    }


    private List<List<List<Double>>> instant_optim(double K_lookahead, double deltaDownload,
                                                   List<Double> p_ij, List<List<List<Double>>> s_ijl, double Ct, List<List<List<Double>>> buf_it,
                                                   List<Double> j_ti, double Bmin, double Bmax){ //give [x_ijl]=
        List<List<List<Double>>> resultx_ijl = s_ijl;

        nb_of_tiles = s_ijl.get(0).get(0).size();
        nb_of_segments = s_ijl.get(0).size();
        nb_of_levels = s_ijl.size();
        List<List<List<Double>>> x_ijl = Matrix.create3DMatrix(nb_of_tiles,nb_of_segments,nb_of_levels,0);


        //-- temporary transform
        List<List<List<Double>>> buf_tmp = buf_it;
        Matrix.setValueInMatrix3DToElementsSmallerThan(buf_tmp,(nb_of_segments+3),1);
        Matrix.setValueInMatrix3DToElementsEqualTo(buf_tmp,(nb_of_segments+3),0);

        //-- end of temporary transform

        List<List<Double>> test = Matrix.sumOfEachMatrix(buf_tmp);
        List<Double> ind_nonfullbuf = Matrix.getValuesSmallerThan(Matrix.vectorSumOfEachRowInMatrix(test), Bmax);
        double nb_of_nonfulltiles = ind_nonfullbuf.size();


        double j_t = Matrix.minVector(j_ti);
        List<Double> j_imin = Matrix.createVector(nb_of_tiles,0);
        List<Double> bufsize_i = Matrix.vectorSumOfEachRowInMatrix(Matrix.sumOfEachMatrix(buf_tmp));

        //-- Initialize with highest qualities on as many next segments as possible
        //with buf size(i)>=Bmax

        for (int ind=0; ind<nb_of_nonfulltiles; ind++){
            int i = ind_nonfullbuf.get(ind).intValue()-1;

            for(int y= j_ti.get(i).intValue()-1;
                y< Math.min(j_t+K_lookahead-1, j_ti.get(i)-1+Bmax-bufsize_i.get(i)) -1;
                y++){
                x_ijl.get(nb_of_levels).get(y).set(i, 1.);
            }

            double newValue = j_ti.get(i)+(Bmin -bufsize_i.get(i))/deltaDownload;
            j_imin.set(i, newValue);
            j_imin.set(i, Math.min(j_t +K_lookahead -1, j_imin.get(i)));
        }


        List<List<List<Double>>>  reqbws = Matrix.multiplyElementByElementMatrix3D(x_ijl,s_ijl);
        double reqbw = Matrix.matrix3Dsum(reqbws);

        List<Double> psorted_i = Matrix.cloneVector(p_ij);
        Collections.sort(psorted_i);

        List<Double> indsorted_i = Matrix.cloneVector(ind_nonfullbuf);


        for (int i=0; i<psorted_i.size(); i++) {
            for (int index=0; index<p_ij.size(); index++) { // get the old index of where the sorted value was
                if(psorted_i.get(index).equals(p_ij.get(i))){
                    indsorted_i.set(index, ind_nonfullbuf.get(index));
                }
            }
        }


        //-- If needed, decrease quality on most future segs then with lowest prob
        //to be watched, still satisfing buf size(i)>=Bmin at the end of dl

        int j = (int) (j_t+K_lookahead-1);
        int l, l_current;
        while (reqbw > Ct*deltaDownload && j>=j_t){
            for (int indsorted=0; indsorted<nb_of_nonfulltiles; indsorted++){
                int i = ind_nonfullbuf.get(indsorted_i.get(indsorted).intValue()).intValue();

                l_current = 0;//TODO : find(x_ijl(i,j,:));

                if (j > j_imin.get(i)){
                    l=l_current-1;
                    x_ijl.get(l_current).get(j).set(i,0.);
                    if (l>=1){
                        x_ijl.get(l).get(j).set(i,1.);
                    }
                } else if (j >= j_ti.get(i)){
                    if (j <= j_imin.get(i) && l_current>1){
                        l=l_current-1;
                        x_ijl.get(l_current).get(j).set(i,0.);
                        x_ijl.get(l).get(j).set(i,1.);
                    }
                }
                reqbws = Matrix.multiplyElementByElementMatrix3D(x_ijl,s_ijl);
                reqbw = Matrix.matrix3Dsum(reqbws);
                if (reqbw <= Ct*deltaDownload)
                    break;
            }
            j = j-1;

            //---Test---//
            if (Matrix.matrix3Dsum(x_ijl)>1)
                System.out.println("more than 1 level per chunk");
        }

        //-- If this is not sufficient, then break constrain buf size(i)>=Bmin and
        //do not dl later segs, in order of p_ij
        if (j<j_t){
            j = (int) (j_t+K_lookahead-1);
            while (j>=j_t && reqbw>Ct*deltaDownload){
                int i;
                for (int indsorted=1; indsorted<nb_of_nonfulltiles; indsorted++){// at least the most important tile downloaded (but for all K_lookahead -> change?)
                    i = ind_nonfullbuf.get(indsorted_i.get(indsorted).intValue()).intValue();

                    //TODO : x_ijl(i,j,:) = zeros(size(x_ijl(i,j,:)));
                    /**
                    for (List<List<Double>> m :
                            x_ijl) {
                        x_ijl.get(i).set(j, x_ijl.get(i).get(j) );
                    }**/

                    reqbws = Matrix.multiplyElementByElementMatrix3D(x_ijl,s_ijl);
                    reqbw = Matrix.matrix3Dsum(reqbws);
                    if (reqbw <= Ct*deltaDownload){
                        break;
                    }
                }
                j=j-1;
            }
        }
        return resultx_ijl;
    }

}