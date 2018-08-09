import java.util.ArrayList;
import java.util.List;

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


        List<List<Double>> bw_trace = Matrix.createMatrix(((int) floor(duration/seg_duration)),1, 1);
        Matrix.multiplyMatrix(cst_bw,bw_trace);

        //---Streaming start---//
        double time_user=0; // counted in number of segments
        double time_video=0; // counted in number of segments played out so far
        int seg_played=0;

        List<List<List<Double>>> buf_it = Matrix.create3DMatrix(Bmax+2,nb_of_tiles,nb_of_levels,1);
        Matrix.multiply3D(-3, buf_it); // stores the segments index, j, 1st level 0 means no seg
        int  nb_of_startupseg = (int) floor(Bmin/2f);

        List<List<Double>> page = Matrix.getPage(buf_it, 1);
        for(int col=0; col < nb_of_startupseg; col++){ //iterate throught the columns from 0 to nb_of_startupseg
            int currentNumber = 1;
            for(int row=0; row < nb_of_startupseg; row++){
                Matrix.setMatrix3DValue(buf_it, row,col,1, currentNumber);
                currentNumber++;
            }
        }
        // Startup period
        List<List<Double>> played_qualities = Matrix.createMatrix(nb_of_tiles,nb_of_segments,0);
        List<List<Double>> dlded_size = Matrix.getUnderMatrixColumnIteration( Matrix.getPage(s_ijl,1), 1 , nb_of_startupseg);
        double dlded_size_sum = Matrix.matrixSum(dlded_size);
        System.out.println(dlded_size_sum);
        time_user = computeDlDelay(dlded_size_sum,bw_trace,0); // startup delay

        /***
        //-- update bw estimate
        double truedl_rate=bw_trace(1:time_user);
        double Ct = estimateNextBtw(1,truedl_rate);

        List<Double> FoV_xy = Matrix.createVector(2,0);
        double total_stalltime=0;

        while (seg_played<nb_of_segments){
            //-- prepare for dl decision
            FoV_xy = generateNewFOV(FoV_xy);
            List<Double> dist_FoVtiles = sum((repmat(FoV_xy,nb_of_tiles,1)-tile_centers_xy).^2,2);
            List<Double> p_ij=(Matrix.maxVector(dist_FoVtiles)-dist_FoVtiles)/sum(Matrix.maxVector(dist_FoVtiles)-dist_FoVtiles);
            List<Double> j_ti=min(buf_it,[],3);
            j_ti=min(j_ti,[],2);
            j_ti(j_ti==(nb_of_segments+3))=time_video;
            j_ti=j_ti+1;
            double j_t = Matrix.minVector(j_ti);

            //-- make dl decision
            List<List<List<Double>>> x_ijl = instant_optim(K_lookahead,deltaDownload,p_ij,s_ijl,Ct,buf_it,j_ti,Bmin,Bmax);

            //-- buffers states after dl finished (right before next download attempt)

            dlded_size = x_ijl.*s_ijl;
            dlded_size_sum = sum(dlded_size(:));
            double time_to_dl = computeDlDelay(dlded_size_sum,bw_trace,time_user); //counted in segments
            //-- temporary transform
            List<List<List<Double>>> buf_tmp = buf_it;
            buf_tmp(buf_tmp < (nb_of_segments+3))=1;
            buf_tmp(buf_tmp == (nb_of_segments+3))=0;
            //-- end of temporary transform
            double cursizebuf_min= Math.min(sum(sum(buf_tmp,3),2));
            double stall_time = Math.max(0,time_to_dl-cursizebuf_min);
            total_stalltime=total_stalltime+stall_time;
            updateBufferStates(buf_it,played_qualities,x_ijl,j_ti,time_to_dl,K_lookahead,time_video,Bmax);

            //-- update bw estimate
            truedl_rate=bw_trace(time_user:time_user+time_to_dl-1);
            Ct = estimateNextBtw(Ct,truedl_rate);

            //-- update times
            time_user = time_user + deltaDownload + Math.ceil(stall_time/seg_duration);
            time_video = time_video+time_to_dl-stall_time;
        }**/

    }


    private void updateBufferStates(double buf_it,double played_qualities,double x_ijl,double j_ti,
                                    double time_to_dl,double K_lookahead,double time_video,double Bmax){ //update buf_it and played_qualities
        /***
        nb_of_tiles=size(x_ijl,1);
        nb_of_segments=size(x_ijl,2);
        j_t=min(j_ti);
        //-- temporary transform
        buf_tmp=buf_it;
        buf_tmp(buf_tmp<(nb_of_segments+3))=1;
        buf_tmp(buf_tmp==(nb_of_segments+3))=0;
        //-- end of temporary transform
        bufsize_i=matrixSum(matrixSum(buf_tmp,3),2);

        //-- first fill: everything that was scheduled for dl
        buf_tmp=buf_it;
        for i=1:nb_of_tiles,
                ind_placeinbuf=1;
        for j=j_ti(i):j_t+K_lookahead-1,
                l=find(x_ijl(i,j,:)==1);
        if ~isempty(l),
                buf_tmp(i,bufsize_i(i)+ind_placeinbuf,l)=j;
        ind_placeinbuf=ind_placeinbuf+1;
        end
                end
        end;

        //-- then drain: everything that has been read in max(deltaDownload,time_to_dl)
        for played_seg=time_video+1:time_video+time_to_dl,
        for i=1:nb_of_tiles,
                ind_l=find(buf_tmp(i,played_seg-time_video,:)==played_seg),
        if length(ind_l)>1 || buf_tmp(i,played_seg-time_video,ind_l)~=played_seg
        error('played_seg index');
        end
        played_qualities(i,played_seg)=ind_l;
        end
                end
        buf_it=buf_tmp(:,time_to_dl+1:min(time_to_dl+Bmax+2,size(buf_tmp,2)),:);
        buf_it(:,size(buf_it,2)+(1:Bmax+2-size(buf_it,2)),:)=nb_of_segments+3;
     ***/
    }

    private List<Double> generateNewFOV(List<Double>[] FoV_xy){
        List<Double> newFov = Matrix.createVector(2,0);
        //TODO : newFov = FoV_xy + range(1,2); // range return an array of 2 randoms numbers
        //TODO: FoV_xy(FoV_xy>1)=FoV_xy(FoV_xy>1)-1;
        return newFov;
    }

    private double estimateNextBtw(double Ct, double truedl_rate) { //update Ct
        double newCt = 0;
        //TODO
        /**
        float alpha = 1;
        for t=1:length(truedl_rate),
                Ct=(1-alpha)*Ct+alpha*truedl_rate(t);
        end
         **/
        return newCt;
    }

    private double computeDlDelay(double dlded_size,List<List<Double>> bw_trace,double start_time_user){ //update time_user
        double delay = 0;
        //TODO
        //size_dldable=cumsum(bw_trace(start_time_user+1:end));
        //time_dlded=find((size_dldable>=dlded_size));
        //time_user=time_dlded(1);
        return delay;
    }


    private List<List<List<Double>>> instant_optim(double K_lookahead, double deltaDownload,
                                                   List<Double> p_ij, List<List<List<Double>>> s_ijl, double Ct, List<List<List<Double>>> buf_it,
                                                   List<Double> j_ti, double Bmin, double Bmax){ //give [x_ijl]=
        List<List<List<Double>>> resultx_ijl = s_ijl;
        /***
        nb_of_tiles=size(s_ijl,1); nb_of_segments=size(s_ijl,2); nb_of_levels=size(s_ijl,3);
        x_ijl=zeros([nb_of_tiles,nb_of_segments,nb_of_levels]);
        //-- temporary transform
        buf_it(buf_it<(nb_of_segments+3))=1;
        buf_it(buf_it==(nb_of_segments+3))=0;
        //-- end of temporary transform
        ind_nonfullbuf=find(matrixSum(matrixSum(buf_it,2),3)<Bmax);
        nb_of_nonfulltiles=length(ind_nonfullbuf);

        j_t=min(j_ti);
        j_imin=zeros(1,nb_of_tiles);
        bufsize_i=matrixSum(matrixSum(buf_it,3),2);

        //-- Initialize with highest qualities on as many next segments as possible
        //with buf size(i)>=Bmax
        for ind=1:nb_of_nonfulltiles,
                i=ind_nonfullbuf(ind);
        x_ijl(i,j_ti(i):min(j_t+K_lookahead-1,j_ti(i)-1+Bmax-bufsize_i(i)),nb_of_levels)=1;
        j_imin(i)=j_ti(i)+(Bmin-bufsize_i(i))/deltaDownload;
        j_imin(i)=min(j_t+K_lookahead-1,j_imin(i));
        end

                reqbw=x_ijl.*s_ijl;
        reqbw=matrixSum(reqbw(:));

        [psorted_i,indsorted_i]=sort(p_ij(ind_nonfullbuf),'ascend');

        //-- If needed, decrease quality on most future segs then with lowest prob
        //to be watched, still satisfing buf size(i)>=Bmin at the end of dl
        j=j_t+K_lookahead-1;
        while reqbw>Ct*deltaDownload && j>=j_t,
        for indsorted=1:nb_of_nonfulltiles,
                i=ind_nonfullbuf(indsorted_i(indsorted));
        l_current=find(x_ijl(i,j,:));
        if j>j_imin(i),
                l=l_current-1;
        x_ijl(i,j,l_current)=0;
        if l>=1,
                x_ijl(i,j,l)=1;
        end
        elseif j>=j_ti(i),
        if j<=j_imin(i) && l_current>1,
                l=l_current-1;
        x_ijl(i,j,l_current)=0;
        x_ijl(i,j,l)=1;
        end
                end
        reqbw=x_ijl.*s_ijl;
        reqbw=matrixSum(reqbw(:));
        if reqbw<=Ct*deltaDownload,
        break;
        end
                end
        j=j-1;
        //---Test---//
        if matrixSum(x_ijl,3)>1,
                error('more than 1 level per chunk');
        end
                end

        //-- If this is not sufficient, then break constrain buf size(i)>=Bmin and
        //do not dl later segs, in order of p_ij
        if (j<j_t){
            j = j_t+K_lookahead-1;
            while (j>=j_t && reqbw>Ct*deltaDownload){
                for (int indsorted=1; indsorted<nb_of_nonfulltiles; indsorted++){// at least the most important tile downloaded (but for all K_lookahead -> change?)
                    i = ind_nonfullbuf(indsorted_i(indsorted));
                    x_ijl(i,j,:) = zeros(size(x_ijl(i,j,:)));
                    reqbw=x_ijl.*s_ijl;
                    reqbw=matrixSum(reqbw(:));
                    if (reqbw <= Ct*deltaDownload){
                        break;
                    }
                }
                j=j-1;
            }
        }
     ***/

        return resultx_ijl;
    }

}