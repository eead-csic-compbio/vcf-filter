/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package vcfreader;

import picard.vcf.SplitVcfs;
import htsjdk.tribble.index.*;
import java.io.IOException;
import java.io.*;

import htsjdk.tribble.index.AbstractIndex;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.tribble.util.LittleEndianOutputStream;


    

   
/**
 *
 * @author Edu
 */
public class Split extends SplitVcfs {
   public String ifile;


    public static void writeTribbleIndex(Index idx, String idxFile) throws IOException {
        LittleEndianOutputStream stream = null;
        try {
            stream = new LittleEndianOutputStream(new BufferedOutputStream(new FileOutputStream(idxFile)));
            idx.write(stream);
        } catch (Exception e) {
           
        // Delete output file as its probably corrupt
            File tmp = new File(idxFile);
            if (tmp.exists()) {
                tmp.delete();
            }
        } finally {
            if (stream != null) {
                stream.close();
            }
        }
    }
    
    public void init() {
    
    int binSize = 16000;

    File inputFile = new File(ifile);
    VCFCodec codec = new VCFCodec();

     String idxFile=ifile+".idx";
     String[] vector=ifile.split("/");
     int ultimo=vector.length-1;
     String[] vector2=vector[ultimo].split("\\.");
     String recomp="";
     for(int i=0;i<vector.length-1; i++){
          recomp=recomp+vector[i]+"/";
     }
     
     recomp=recomp+""+vector2[0];
 

  
    AbstractIndex idx= IndexFactory.createLinearIndex(inputFile, codec, binSize); 
   
      try{
         writeTribbleIndex(idx,idxFile); 
         File fileo1 = new File(recomp+".snp.vcf");
         File fileo2 = new File(recomp+".indel.vcf");

        this.INPUT=inputFile;
        
        this.SNP_OUTPUT=fileo1;
        this.INDEL_OUTPUT=fileo2;
        this.STRICT=false;

        super.doWork();
         } catch (Exception e){
            
         } 
       
      
    }
    
    public static void main(String[] args) {
         
     new Split().init();
     
    }
  
        
    
}

