/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package vcfreader;

/**
 *
 * @author Edu
 */

import htsjdk.tribble.index.AbstractIndex;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.util.LittleEndianOutputStream;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.vcf.VCFFileReader;
import java.io.File;
import java.io.IOException;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterFactory;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Scanner;
import java.util.Set;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.ParseException;

public class VcfReader {

    public int num;
    public String ifile;
    public String pathin="";
    public String pathout="";

    public VCFFileReader VCFreader;
    public Split split;
    public int DPG;
    public int DPS;
    public String variante;
    public String sample;
    public String crom;
    public int position;
    public int posfirstint;
    public int possecondint;
    public double minHet;
    public int minGQ;
    public int Size;
    public Set set;
    public static CommandLine cmd;
    
    public static org.apache.commons.cli.Options options;
    
    public int total;
    public double frecrar;
    public double NData;
    public VariantContextWriter vcfwriter;
    public EnumSet<Options> DEFAULT_OPTIONS=
            EnumSet.of(Options.INDEX_ON_THE_FLY, 
            Options.ALLOW_MISSING_FIELDS_IN_HEADER, 
            Options.WRITE_FULL_FORMAT_FIELD);
    
    
    
    
    public void MinDPGen(){
               
        Iterator<VariantContext> iter=VCFreader.iterator();
        while(iter.hasNext()){
        VariantContext variant= iter.next();
        int DPV=variant.getAttributeAsInt("DP", 0);
        if(DPV>DPG){
            variante=variante+variant.toString()+"\n";
            vcfwriter.add(variant);
        }
        }
       
    }
    public void MinDPSample(){
 
       total=VCFreader.getFileHeader().getNGenotypeSamples();
        Iterator<VariantContext> iter=VCFreader.iterator();
        while(iter.hasNext()){
            VariantContext variant= iter.next();
            boolean minimo=false;
            for(int i=0; i<total;i++){
            int DP= variant.getGenotype(i).getDP();
               if(DPS>DP){
                   minimo=true;
               } 
            }
            if(minimo==false){
              variante=variante+variant.toString()+"\n";
              vcfwriter.add(variant);
            }
        }
     
     }   
     


     public void MinHet(){
         
         total=VCFreader.getFileHeader().getNGenotypeSamples();
         double maxhet=(double)minHet*(double)total/100;
          
         VCFreader.iterator().forEachRemaining(variantcontext -> {
             int numhet=0;
             
             for (int i=0; i<total; i++){
                 if(variantcontext.getGenotype(i).isHet()){
                 numhet=numhet+1;
                 }                
             }
             if (maxhet>numhet){
                 vcfwriter.add(variantcontext);
             }      
             
         });

         
        }
     
    public void Variantrare(){
           
           total=VCFreader.getFileHeader().getNGenotypeSamples();
            HashMap h=new HashMap();
           VCFreader.iterator().forEachRemaining(variantcontext -> {
            double totalalleles=0;
            h.clear();
               for (int i=0; i<total; i++){
                   for (Allele o: variantcontext.getGenotype(i).getAlleles()){
                   totalalleles=totalalleles+1;    
                   String b=o.getBaseString();
                   if(h.containsKey(b)){
                   h.put(b, (int)h.get(b)+1 );
                   }
                   else {
                   h.put(b, 1);
                   }
                   
                   
                   }

                 }
               boolean pass=true;
               Iterator iter=h.keySet().iterator();
               while(iter.hasNext()){        
               String c=(String)iter.next();
               double frec=(int)h.get(c)/totalalleles;
               if(frec<frecrar){
                   pass=false;
               }
               }
               if(pass==true){
               vcfwriter.add(variantcontext);
               }
             

             }); 
         
        } 

    public void MissingData(){
        
       
       total=VCFreader.getFileHeader().getNGenotypeSamples();
       NData=total*NData/100;


       VCFreader.iterator().forEachRemaining(variantcontext -> {
       int num=0;
       for(int i=0; i<total;i++){
      
        int DP2= variantcontext.getGenotype(i).getDP();
           if(DP2==0){
               num=num+1;
           } 
           
        }           
        if(num<NData){
            vcfwriter.add(variantcontext);
        }
        });
       
       
       }
    
    public void NumBiallelic(){
           
           VCFreader.iterator().forEachRemaining(variantcontext -> {
    
            if(variantcontext.isSNP() && variantcontext.isBiallelic()){
                vcfwriter.add(variantcontext);
            
           } 
        });
       }
    
    public void CreateVCF()throws IOException {
        VCFreader =new VCFFileReader(new File(pathin));
        if(pathout.isEmpty()){
                vcfwriter=VariantContextWriterFactory.createVcf(null, System.out,
                VCFreader.getFileHeader().getSequenceDictionary(),DEFAULT_OPTIONS);
        }
        else {
                FileOutputStream outputstream= new FileOutputStream(new File(pathout));
                vcfwriter=VariantContextWriterFactory.createVcf(new File(pathout), outputstream,
                VCFreader.getFileHeader().getSequenceDictionary(),DEFAULT_OPTIONS);

        }



        vcfwriter.writeHeader(VCFreader.getFileHeader());

        
    }
    
    public void FindSamVar(){
       
       Iterator<VariantContext> iter=VCFreader.iterator();
       boolean find=false;

        while(iter.hasNext()){
            VariantContext variant= iter.next();
            int num=variant.getEnd();
            String name=variant.getContig();
            if(num==position && name.equals(crom)){
                variante=variant.getGenotype(sample).toString();
                find=true;
                break;
             
            }    
        }
        if(find==false){
        System.out.println("your sample doesn't exit");
        }

    }
          
    
    
    public void SelectGenotype() throws IOException {
       
         VCFHeader header= new VCFHeader(VCFreader.getFileHeader().getMetaDataInInputOrder(), set);
         FileOutputStream outputstream= new FileOutputStream(new File(pathout));
         vcfwriter=VariantContextWriterFactory.createVcf(new File(pathout), 
                outputstream, VCFreader.getFileHeader().getSequenceDictionary(),DEFAULT_OPTIONS);
     
          vcfwriter.writeHeader(header); 
            
        Iterator<VariantContext> iter= VCFreader.iterator();
        while(iter.hasNext()){
        VariantContext variant=iter.next();
         VariantContext vc2=variant.subContextFromSamples(set);
                vcfwriter.add(vc2);
        }
    
    }

    public void VariantInter(){

       
       Iterator<VariantContext> iter=VCFreader.iterator();


       while(iter.hasNext()){
            VariantContext variant= iter.next();
             String id= variant.getContig();
    
    
         if (variant.isSNP() && variant.getEnd()>posfirstint && variant.getEnd()<possecondint && id.equals(crom)){
             vcfwriter.add(variant);
             variante=variante+variant.toString()+"\n";

          }
         if(variant.getEnd()>possecondint && id.equals(crom)){
         break;
         }
        }

    }
    
    public void Split(){
        Split splits= new Split();
        splits.ifile=pathin;
        splits.init();
    }
    
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
    
    
    public void createidx() {
    
    int binSize = 16000;
    File inputFile = new File(pathin);
    VCFCodec codec = new VCFCodec();

    String idxFile=pathin+".idx";

  
    AbstractIndex idx= IndexFactory.createLinearIndex(inputFile, codec, binSize); 
   
     try{
         writeTribbleIndex(idx,idxFile);      
    }catch (Exception e){}
     
    }
    
    public void BestQUALinKb(){

        Iterator<VariantContext> ita= VCFreader.iterator();  
        ita.next();
        VariantContext variantref= ita.next();

        int posref= variantref.getEnd()+Size;
        String cromo=variantref.getContig();
        double qual=variantref.getPhredScaledQual();
        int posvariant=0;
        int i=0;
        VariantContext max=variantref;
        while(ita.hasNext()){
            
            VariantContext variant =ita.next();
            posvariant=variant.getEnd();
            
                if(cromo.equals(variant.getContig())==false){
                    vcfwriter.add(max);
                    cromo=variant.getContig();
                    posref=variant.getEnd()+Size;
                    qual=variant.getPhredScaledQual();
                    max=variant;
                
                }
            
                if (posvariant<posref && variant.getContig().equals(cromo) && variant.getPhredScaledQual()>qual){
                    qual=variant.getPhredScaledQual();
                    max=variant;
                    
                }

                if(posvariant>posref && variant.getContig().equals(cromo)){
                    vcfwriter.add(max);
                    posref=variant.getEnd()+Size;
                    qual=variant.getPhredScaledQual();
                    max=variant;
                
                }

             i++;

        }
    
    }
    
    public void MinGQkb(){

        total=VCFreader.getFileHeader().getNGenotypeSamples();
        ArrayList<Integer> array= new ArrayList();

        Iterator<VariantContext> ita= VCFreader.iterator();  
        ita.next();
        VariantContext variantref= ita.next();

        int posref= variantref.getEnd()+Size*1000;
        String cromo=variantref.getContig();
        int num=0;
        int posvariant=0;
        while(ita.hasNext()){
            
            VariantContext variant =ita.next();
            boolean find=false;
            posvariant=variant.getEnd();
                if (posvariant<posref && variant.getContig()==cromo){
                    for (int i=0; i<total;i++){
                        if(variant.getGenotype(i).getGQ()<minGQ){
                            find=true;      
                        }
                        
                    } 
                    if(find==false){
                    num+=1;
                    }
                              
                }
                else {
                posref=variant.getEnd()+Size*1000;
                cromo=variant.getContig();
                array.add(num);
                num=0;
                
                }
        }
        variante=""+array.size();
  
    
    }
    
    public void menu() throws IOException {
        
        Scanner sc= new Scanner(System.in);   
        System.out.println("First, type the path to your INPUT file");             
        pathin=sc.next();
        System.out.println("Now, type the path to your OUTPUT file");
        pathout=sc.next();
        int opt;
        do {
        System.out.println("-----------------------------");
        System.out.println("Options Menu in VCF File");
        System.out.println("-----------------------------");
        System.out.println("Option 1- Filter by minimun General DP in VariantContext");
        System.out.println("Option 2- Filter by minimun  DP in each Genotype");
        System.out.println("Option 3- Filter by maximun of missing data of genotypes");
        System.out.println("Option 4- Get a sample in a specific position");
        System.out.println("Option 5- Obtain number of SNPS with minimun  sample GQ in each Genotype in a specific distance in kb");
        System.out.println("Option 6- Select SNPs in an interval in a concrete cromosome");
        System.out.println("Option 7- Filter by Biallelic SNPs");
        System.out.println("Option 8- Filter by minimum allele frequency");
        System.out.println("Option 9- Filter by % maximun of Heterozigotes samples");
        System.out.println("Option 10- Split VCF File in to SNPs and Indels");
        System.out.println("Option 11- Select samples and generate new vcf file");
        System.out.println("Option 12- Select the best quality variant in your desired distance");
        System.out.println("Press 0 to leave");
        System.out.println("----------------------------------------------------");   
        System.out.println("Choose the option number");
        opt=sc.nextInt();
        
        switch (opt) {

            case 1:
                System.out.println("Type the minimum General DP");
                DPG=sc.nextInt();
                CreateVCF();
                MinDPGen();

            break;
            
            case 2:
                System.out.println("Type the minimum sample DP");
                DPS=sc.nextInt();
                CreateVCF();
                MinDPSample();
            
            break;
                
            case 3:
                System.out.println("Type the % of maximum missing data per variantcontext");
                NData=sc.nextInt();
                CreateVCF();
                MissingData();
            
            break;
            
            case 4:
                System.out.println("Type the sample to find");
                sample=sc.next();
                System.out.println("Type the specific position");
                position=sc.nextInt();
                System.out.println("Type the specific chromosome");
                crom=sc.next();
                FindSamVar();
                System.out.println(variante);
            
            break;
                        
            case 5:
                System.out.println("Type the minimun sample GQ");
                minGQ=sc.nextInt();
                System.out.println("Type the distance in Kb");
                Size=sc.nextInt();
                CreateVCF();
                MinGQkb();
                System.out.println(variante);
            
            break;
            
            case 6:
                System.out.println("Introduce the first position");
                posfirstint=sc.nextInt();
                System.out.println("Introduce the last position");
                possecondint=sc.nextInt();
                System.out.println("Introduce the specific cromosome");
                crom=sc.next();
                CreateVCF();
                VariantInter();
            
            break;
            
            case 7:
                CreateVCF();
                NumBiallelic();       
            
            break;
            
            case 8:
                System.out.println("Introduce the number of minimum allele frequency. Use , instead of .");
                frecrar=sc.nextDouble();
                CreateVCF();
                Variantrare();
            
            break;
            
            case 9:
                System.out.println("Introduce % of maximum heterozigotes samples");
                minHet=sc.nextDouble();
                CreateVCF();
                MinHet();
            
            break;
            
            case 10:
                Split();
              
            break;
            
            case 11:
                System.out.println("Write selected samples separated by ;");
                String samples1=sc.next();
                set= new HashSet();
                String [] samples2=samples1.split(";");
                System.out.println(samples2.length);
                for (int i=0; i<samples2.length; i++){
                    System.out.println(samples2[i]);
                set.add(samples2[i]);
                }
                SelectGenotype();
            break;
            
            case 12:
                System.out.println("Write your desired distance");
                Size=sc.nextInt();
                CreateVCF();
                BestQUALinKb();
            break;
            
            case 0:
            break;    
            
            default:
                System.out.println("Introduce correct option");    
            break;
   
        }         
        } while (opt!=0);
    
    
    
    }
    
    public void ReadOptions() throws IOException {
       


          
        if (cmd.hasOption("input")){
            pathin=cmd.getOptionValues("input")[0];
            createidx();
        }
        
        if (pathin.isEmpty()){  
        }
        else {
        VCFreader =new VCFFileReader(new File(pathin));
        }
           

        if (cmd.hasOption("output")){
          pathout=cmd.getOptionValues("output")[0];
        } 
        
        if(cmd.hasOption("DP")) {  
                DPG=Integer.parseInt(cmd.getOptionValues("DP")[0]);
                CreateVCF();
                MinDPGen();
                System.out.println("hola");
                
            }
        if (cmd.hasOption("sDP")) {     
                DPS=Integer.parseInt(cmd.getOptionValues("sDP")[0]);
                CreateVCF();
                MinDPSample(); 
        }
        
        if (cmd.hasOption("missing")) {     
                NData=Integer.parseInt(cmd.getOptionValues("missing")[0]);
                CreateVCF();
                MissingData();
        }
        if(cmd.hasOption("crom")){
           crom=cmd.getOptionValues("crom")[0];
        }
        
        if (cmd.hasOption("interval")){
                posfirstint=Integer.parseInt(cmd.getOptionValues("interval")[0]);
                possecondint=Integer.parseInt(cmd.getOptionValues("interval")[1]);
                CreateVCF();
                VariantInter();
        }
        
        if(cmd.hasOption("pos")){
            position=Integer.parseInt(cmd.getOptionValues("pos")[0]);
        }
        
        if (cmd.hasOption("sample")){
                        
                set= new HashSet();
                String [] samples2=cmd.getOptionValues("sample")[0].split(",");
                System.out.println(samples2.length);
                for (int i=0; i<samples2.length; i++){
                    System.out.println(samples2[i]);
                set.add(samples2[i]);
                }
                SelectGenotype();
        
        }
        
        if (cmd.hasOption("call")){
                sample=cmd.getOptionValues("call")[0];
                FindSamVar();
                System.out.println(variante);
        
        }
        if (cmd.hasOption("bi")){
                CreateVCF();
                NumBiallelic(); 
        
        }
        
        if (cmd.hasOption("MAF")){
                frecrar=Double.parseDouble(cmd.getOptionValues("MAF")[0]);
                CreateVCF();
                Variantrare();
        
        }
        
        if (cmd.hasOption("maxHet")){
                minHet=Double.parseDouble(cmd.getOptionValues("maxHet")[0]);
                CreateVCF();
                MinHet();
        }
        
        if (cmd.hasOption("nr")){
               Size=Integer.parseInt(cmd.getOptionValues("nr")[0]);
                CreateVCF();
                BestQUALinKb();
        
        }
        
        if(cmd.hasOption("menu")){
            menu();
        
        }
        
        if (cmd.hasOption("split")){
            Split();
        
        }
        

        
    }
    
    public void CreateOptions(){
        
        
        Option input = OptionBuilder.withArgName( "path to input uncompressed VCF file" )                    
                                .hasArgs(1)
                                .withDescription(  "optional, by default reads from STDIN" )
                                .create( "input" );
        
        Option output = OptionBuilder.withArgName( "path to output VCF file" )                    
                                .hasArgs(1)
                                .withDescription("optional, by default prints to STDOUT")
                                .create( "output" );        
        
        Option DP = OptionBuilder.withArgName( "integer" )
                                .hasArgs(1)
                                .withDescription("Minimum overall DP of each variant (row) in the input VCF")
                                .create( "DP" );
        
        Option sDP = OptionBuilder.withArgName( "integer" )
                                .hasArgs(1)
                                .withDescription("Minimum DP of each sample")
                                .create( "sDP" );
        
        Option missing = OptionBuilder.withArgName( "double" )
                                .hasArgs(1)
                                .withDescription("Allowed missing data in each variant percentage \n example: -missing 5 ")
                                .create( "missing" );
        
        Option interval = OptionBuilder.withArgName( "integer" )
                                .hasArgs(2)
                                .withDescription("Select a particular interval requires -crom \n example: -interval 1000 1000 ")
                                .create( "interval" );
        
        Option crom = OptionBuilder.withArgName( "String" )
                                .hasArgs(1)
                                .withDescription( "Select the cromosome of your interest" )
                                .create( "crom" );
        
        Option pos = OptionBuilder.withArgName( "integer" )
                                .hasArgs(1)
                                .withDescription("Select the position of your interest")
                                .create( "pos" );
        
        Option sample = OptionBuilder.withArgName( "String: samples names" )
                                .hasArgs()
                                .withDescription("Select a subset of samples\n example: -sample sample1,sample2" )
                                .create( "sample" );
        
        Option call = OptionBuilder.withArgName( "String: sample name" )
                                .hasArgs(1)
                                .withDescription("requires -interval, -sample \n example: -sample sample1 -crom Bd1 -pos 1000" )
                                .create( "call" );
        
        Option bi = new Option( "bi", "Select only biallelic genotypes " );
        
        Option MAF = OptionBuilder.withArgName("double")
                                .hasArgs(1)
                                .withDescription(  "Minimum allele frequency [0-1]\n example: -MAF 0.05" )
                                .create( "MAF" );
        
        Option maxHet = OptionBuilder.withArgName("double")
                                .hasArgs(1)
                                .withDescription("Maximum Heterozigous samples\n percentage, example: -maxHet 90" )
                                .create( "maxHet" );
        
        Option nr = OptionBuilder.withArgName("integer in Kb")
                                .hasArgs(1)
                                .withDescription("Select top quality variant in window\n example: -nr 100" )
                                .create( "nr" );

        
        Option split = new Option( "split", "Split in two vcf files, SNPs and Indels" );
        
   
        Option menu = new Option( "menu", "Select helped menu" );

        options.addOption(DP);
        options.addOption(input);
        options.addOption(output);
        options.addOption(sDP);
        options.addOption(missing);
        options.addOption(interval);
        options.addOption(crom);
        options.addOption(pos);
        options.addOption(sample);
        options.addOption(bi);
        options.addOption(MAF);
        options.addOption(maxHet);
        options.addOption(nr);
        options.addOption(call);
        options.addOption(menu);
        options.addOption(split);

        

        
                

    
    }
    
    public void command() throws IOException{
        Scanner sc= new Scanner(System.in);
        System.out.println("This is an application for get specific information in a VCF file");
        System.out.println("In most cases generates a new VCF file in the desirable output rout");

        System.out.println("\033[0;1mOptions\033[0m"+"\n");
        System.out.println("GeneralDP : Introduce the Minimum DP in each variantcontext");
        System.out.println("SampleDP : Introduce the Minimum DP per sample");
        System.out.println("MissingData : Introduce the Maximun Missing data in each variantcontext");
        System.out.println("Usage example: GeneralDP [minDP] [input=your interest VCF file rout] [output=your interest output rout]");
        System.out.println("--------------------------------------------------------------------------------------------------------");
        System.out.println("Find : Introduce the sample in a specific cromosome in specific localitation");
        System.out.println("Usage example: Find [namesample] [namechromosome] [position] [input=your interest VCF file rout]");
        System.out.println("--------------------------------------------------------------------------------------------------------");
        System.out.println("SelectSNPs : Select SNPs in a concrete interval");
        System.out.println("Usage example: SelectSNPs [inters=1000-2000] [namecromosome] [input=your interest VCF file rout] [output=your interest output rout]");
        System.out.println("--------------------------------------------------------------------------------------------------------");       
        System.out.println("Biallelic: Select only the SNPs that are Biallelic");
        System.out.println("Usage example: Biallelic [input=your interest VCF file rout] [output=your interest output rout]");
        System.out.println("--------------------------------------------------------------------------------------------------------");      
        System.out.println("MAF : Select your minimun allele frequency ");
        System.out.println("MaxHet : Select a threshold of maximum Heterozigotes samples");
        System.out.println("Usage example: MAF [MAF] [input=your interest VCF file rout] [output=your interest output rout]");
        System.out.println("--------------------------------------------------------------------------------------------------------");      
        System.out.println("Split: Split your VCF file in SNPs and indels");
        System.out.println("Usage example: Split [input=your interest VCF file rout]");
        System.out.println("--------------------------------------------------------------------------------------------------------");      
        System.out.println("SelectSamples: Generates new vcf file with the sample of your interest, the samples have to separete by ; ");
        System.out.println("Usage example: SelectSamples [samplename1;samplename2;..;] [input=your interest VCF file rout] [output=your interest output rout]");
        System.out.println("--------------------------------------------------------------------------------------------------------");
        System.out.println("BestQuality: Generates new vcf file with the top quality variant in desired distance ");
        System.out.println("Usage example: BestQuality [distance in Kb] [input=your interest VCF file rout] [output=your interest output rout]");
        System.out.println("--------------------------------------------------------------------------------------------------------");
        
       
        String comands=sc.nextLine();
        String [] option=comands.split(" ");
        String input="";
        String output="";
        for (int i=0; i<option.length; i++){
        if (option[i].contains("input")){ 
            String prue[]=option[i].split("=");
            input=prue[1];
                }
        if (option[i].contains("output")){ 
            String prue[]=option[i].split("=");
            output=prue[1];
                }        
        }
        
        if (input.isEmpty()==false){
         pathin=input;
        }
        
      
        if (output.isEmpty()){
        pathout="";
        }
        else {
        pathout=output;
        }

       
         
        
        switch (option[0]) {

            case "GeneralDP":
                
                DPG=Integer.parseInt(option[1]);
                CreateVCF();
                MinDPGen();

            break;
            
            case "SampleDP":
                DPS=Integer.parseInt(option[1]);
                CreateVCF();
                MinDPSample();
            break;
            
            case "MissData":
                NData=Integer.parseInt(option[1]);
                CreateVCF();
                MissingData();
            break;
            
            case "Find":
                sample=option[1];
                crom=option[2];
                position=Integer.parseInt(option[3]);
                FindSamVar();
                System.out.println(variante);
            break;
            
            case "SelectSNPs":
                
                String inter[]=option[1].split("=");
                String inters[]=inter[1].split("-");

                posfirstint=Integer.parseInt(inters[0]);
                possecondint=Integer.parseInt(inters[1]);
                crom=option[2];
                CreateVCF();
                VariantInter();
            break;
            
            case "Biallelic":
                CreateVCF();
                NumBiallelic(); 
            break;
            
            case "MAF":
                frecrar=Double.parseDouble(option[1]);
                CreateVCF();
                Variantrare();
            break;
            
            case "MaxHet":
                minHet=Double.parseDouble(option[1]);
                CreateVCF();
                MinHet();
            break;
            
            case "Split":
                Split();
            break;
            
            case "SelectSamples":
                set= new HashSet();
                String [] samples2=option[1].split(";");
                System.out.println(samples2.length);
                for (int i=0; i<samples2.length; i++){
                    System.out.println(samples2[i]);
                set.add(samples2[i]);
                }
                SelectGenotype();
            break;
            
            case "BestQuality":
                Size=Integer.parseInt(option[1]);
                CreateVCF();
                BestQUALinKb();
            break;
            
            default:
                System.out.println("Introduce one of the execute options or put it in the apropiate format");
            break;    
                
        }

    
    }
    public static void main(String[] args)throws IOException, ParseException {

     
        options = new org.apache.commons.cli.Options();
        
        VcfReader vcffile=new VcfReader();
        
        vcffile.CreateOptions();
        
        CommandLineParser parser = new DefaultParser(); 

        cmd = parser.parse( options, args);
        
        


        vcffile.ReadOptions();

 
        
        HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp("help menu", options);
        
        System.out.println("This utility builds on HTSJDK classes, apache commons CLI classes and supports VCF versions supported therein.");
        System.out.println("Eduardo Candeal, Carlos P Cantalapiedra, Bruno Contreras-Moreira\n" +
"EEAD-CSIC 2016");
    }
    
    
    
}
