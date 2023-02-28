package masterPATH;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Paths;
import java.util.Map;
//import org.apache.commons.cli;
//import org.apache.derby;
import org.apache.commons.cli.*;

/**
 * Main class
 *
 * @author Natalia Rubanova
 */
public class masterPATH {

    /**
     * main method
     *
     * @param args
     * @throws java.io.FileNotFoundException
     * @throws java.lang.InterruptedException
     */
    public static void main(String[] args) throws FileNotFoundException, IOException, InterruptedException {

        String userDirectory = Paths.get("")
                .toAbsolutePath()
                .toString();

        //System.out.println(userDirectory);
//        String hugo_path = userDirectory + File.separator + "files" + File.separator + "hugo_all_17Jun.txt" ;
//        String hugo_entrez_path = userDirectory + File.separator + "files" + File.separator + "HUGO_with_entrezid.txt" ;
//        String entrz_wo_hugo= userDirectory + File.separator + "files" + File.separator + "entrez_wo_hugo.txt" ;
//        String hprd_bioDB = userDirectory + File.separator + "files" + File.separator + "missed_hprd_to_hugo_bioDB.txt" ;
//        String tfacts_crosstable = userDirectory + File.separator + "files" + File.separator + "tfacts_crosstable.txt" ;
      //  String hugo_path = "/files" + File.separator + "hugo_all_17Jun.txt";
      //  String hugo_entrez_path = "/files" + File.separator + "HUGO_with_entrezid.txt";
      //  String entrz_wo_hugo = "/files" + File.separator + "entrez_wo_hugo.txt";
      //  String hprd_bioDB = "/files" + File.separator + "missed_hprd_to_hugo_bioDB.txt";
      //  String tfacts_crosstable = "/files" + File.separator + "tfacts_crosstable.txt";
        
        String hugo_path = "IDconversion" + File.separator + "hugo_all_17Jun.txt";
        String hugo_entrez_path = "IDconversion" + File.separator + "HUGO_with_entrezid.txt";
        String entrz_wo_hugo = "IDconversion" + File.separator + "entrez_wo_hugo.txt";
        String hprd_bioDB = "IDconversion" + File.separator + "missed_hprd_to_hugo_bioDB.txt";
        String tfacts_crosstable = "IDconversion" + File.separator + "tfacts_crosstable.txt";


        DBManager dbutils = new DBManager(hugo_path, hugo_entrez_path, entrz_wo_hugo, hprd_bioDB, tfacts_crosstable);
        NetworkManager nutils = new NetworkManager();
        PathwayManager putils = new PathwayManager();
        CentralityManager outils = new CentralityManager();
        NetworkTopology utils = new NetworkTopology();
        Wrapper wrap = new Wrapper(hugo_path, hugo_entrez_path, entrz_wo_hugo, hprd_bioDB, tfacts_crosstable);
        RandomManager rand = new RandomManager();

        Options options = new Options();

        Option network = new Option("n", "network", true, "Network type: ppi or integrated");
        network.setRequired(true);
        options.addOption(network);

        Option hitlist = new Option("hl", "hitlist", true, "Hit list file name");
        hitlist.setRequired(true);
        options.addOption(hitlist);

        Option finalimpls = new Option("fp", "finalimplementers", true, "Final implementers file name");
        finalimpls.setRequired(true);
        options.addOption(finalimpls);

        Option modeparam = new Option("m", "mode", true, "Mode : paths or network");
        modeparam.setRequired(true);
        options.addOption(modeparam);

        Option pathlistparam = new Option("l", "pathlist", true, "File name for list of pathways ids to create a network file for Cytoscape. One ID on a line");
        pathlistparam.setRequired(false);
        options.addOption(pathlistparam);

        Option maxlength = new Option("ml", "maxlength", true, "Max length of paths for Breadth first algorithm");
        maxlength.setRequired(false);
        options.addOption(maxlength);

        Option minpathlength = new Option("minl", "minpathlength", true, "Min length of paths to report");
        minpathlength.setRequired(false);
        options.addOption(minpathlength);

        Option maxpathlength = new Option("maxl", "maxpathlength", true, "Max length of paths to report");
        maxpathlength.setRequired(false);
        options.addOption(maxpathlength);

        Option mincentrality = new Option("minc", "mincentrality", true, "Min cintrality to test");
        mincentrality.setRequired(false);
        options.addOption(mincentrality);

        Option maxthreads = new Option("t", "threads", true, "Number of threads");
        maxthreads.setRequired(false);
        options.addOption(maxthreads);

        Option randomint = new Option("r", "randominteractions", true, "Number of random interactions");
        randomint.setRequired(false);
        options.addOption(randomint);

        Option pathprefix = new Option("p", "prefix", true, "Prefix for paths' ids");
        pathprefix.setRequired(false);
        options.addOption(pathprefix);

        Option outputfolder = new Option("o", "output", true, "Output folder. Should be specified with -v parameter for docker container ");
        outputfolder.setRequired(true);
        options.addOption(outputfolder);

        CommandLineParser parser = new DefaultParser();
        HelpFormatter formatter = new HelpFormatter();
        CommandLine cmd = null;

        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.out.println(e.getMessage());
            formatter.printHelp("masterPath", options);

            System.exit(1);
        }

        String networkType = cmd.getOptionValue("network");
        String hitList = cmd.getOptionValue("hitlist");
        String finalImplementers = cmd.getOptionValue("finalimplementers");
        String maxLengthString = cmd.getOptionValue("maxlength");
        String minPathLengthString = cmd.getOptionValue("minpathlength");
        String maxPathLengthString = cmd.getOptionValue("maxpathlength");
        String minCentralityString = cmd.getOptionValue("mincentrality");
        String threadsString = cmd.getOptionValue("threads");
        String maxRandomString = cmd.getOptionValue("randominteractions");
        String prefix = cmd.getOptionValue("prefix");
        String outputFolder = cmd.getOptionValue("output");
        String mode = cmd.getOptionValue("mode");
        String pathlist = cmd.getOptionValue("pathlist");
        
        
        
        System.out.println(finalImplementers);
        
        /**************************************************************/
        //Angela: assign value to parameters mannually
        maxLengthString = "100";
        minPathLengthString = "20";
        maxPathLengthString = "80";
        minCentralityString = "3";
        threadsString = "1";
        maxRandomString = "3";
        prefix = "testrunning";
        
        networkType = "integrated";
        outputFolder = "Angela_output";
        mode = "paths";
        finalImplementers = "Examples/Muscle Differentiation/Final implementers";
        hitList = "Examples/Muscle Differentiation/miRNA_hit_list";
        
        /********the above can also be implemented in arguments as follows*******/
        /*
        -n integrated  -hl "/home/bigangela03/eclipse-test/masterPATH/Examples/Muscle Differentiation/miRNA_hit_list"
        -fp "Examples/Muscle Differentiation/Final implementers"
        -o "Angela_output"
        -m paths
        */
        
        /**************************************************************/
        
        
        
        

        // System.out.println(hitList);
        if (!"paths".equals(mode) && !"network".equals(mode)) {
            System.err.println(mode.equals("paths"));
            System.err.println("Mode should be either paths or network");
            System.exit(1);
        }

        if (!"test".equals(hitList)) {
            File hl = new File(hitList);
            if (!(hl.exists() && hl.isFile())) {
                System.err.println("Hit list file does not exist");
                System.exit(1);
            }
        } else {
            System.out.println("Loading test hit list ...  ");
            hitList = "/files" + File.separator + "Hit_genes_test";

        }

        if (!"test".equals(finalImplementers)) {
            File fi = new File(finalImplementers);
            if (!(fi.exists() && fi.isFile())) {
                System.err.println("Final implementer file does not exist");
                System.exit(1);
            }
        } else {
            System.out.println("Loading test final implementers ...  ");
            finalImplementers = "/files" + File.separator + "Final_implementers_test";
        }

        File of = new File(outputFolder);
        if (!(of.exists() && of.isDirectory())) {
            System.err.println("Output folder does not exist");
            System.out.println(outputFolder);
            System.exit(1);
        } else {
            outputFolder = outputFolder + File.separator;
        }

        Network nw = new Network();
        Wrapper wr = new Wrapper(hugo_path, hugo_entrez_path, entrz_wo_hugo, hprd_bioDB, tfacts_crosstable);

        String nwType = "ppi";
        switch (networkType) {
            case "integrated":
//                nw = wr.load_network(
//                        userDirectory + File.separator + "files" + File.separator + "Nodes.updatedApril2019", 
//                        userDirectory + File.separator + "files" + File.separator + "Interactions.updatedApril2019"
//                );
            	/*
                nw = wr.load_network(
                        "/files" + File.separator + "Nodes.updatedApril2019",
                        "/files" + File.separator + "Interactions.updatedApril2019"
                );
                */
            	
            	 nw = wr.load_network(
                         "Networks" + File.separator + "Nodes.updatedApril2019",
                         "Networks" + File.separator + "Interactions.updatedApril2019"
                 );

                nwType = "integrated";
                break;

            case "ppi":
//                nw = wr.load_network(
//                        userDirectory + File.separator + "files" + File.separator + "nodes_hippie_h",
//                        userDirectory + File.separator + "files" + File.separator + "int_hippie_h"
//                );
            	
            	//Angela: nodes_hippie_h can't be found in Examples folder!
                nw = wr.load_network(
                        "/files" + File.separator + "nodes_hippie_h",
                        "/files" + File.separator + "int_hippie_h"
                                            
                );
                nwType = "ppi";
                break;
            default:
                System.err.println("Invalid network type. Please choose between integrated and ppi");
                System.exit(1);
                break;
        }

//        nutils.load_hitlist_and_finalimpl(nw, hitList, finalImplementers, dbutils.hugo);
        int maxLength = 0;
        int minPathLength = 0;
        int maxPathLength = 0;
        int threads = 0;
        int maxRandom = 0;
        int minCentrality = 0;

        if ("paths".equals(mode)) {
            try {
                if ("null".equals(maxLengthString)) {
                    maxLengthString = "na";
                }
                maxLength = Integer.parseInt(maxLengthString);
            } catch (NumberFormatException e) {
                System.err.println("Max length of paths for Breadth first algoritm should be an integer value");
                System.exit(1);
            }

            try {
                if ("null".equals(minPathLengthString)) {
                    minPathLengthString = "na";
                }
                minPathLength = Integer.parseInt(minPathLengthString);
            } catch (NumberFormatException e) {
                System.err.println("Min length of paths for Breadth first algoritm should be an integer value");
                System.exit(1);
            }

            try {
                if ("null".equals(maxPathLengthString)) {
                    maxPathLengthString = "na";
                }
                maxPathLength = Integer.parseInt(maxPathLengthString);
            } catch (NumberFormatException e) {
                System.err.println("Min length of paths to report should be an integer value");
                System.exit(1);
            }

            try {
                if ("null".equals(minCentralityString)) {
                    minCentralityString = "na";
                }
                minCentrality = Integer.parseInt(minCentralityString);
            } catch (NumberFormatException e) {
                System.err.println("Min centrality of paths to report should be an integer value");
                System.exit(1);
            }

            try {
                if ("null".equals(threadsString)) {
                    threadsString = "na";
                }
                threads = Integer.parseInt(threadsString);
            } catch (NumberFormatException e) {
                System.err.println("Threads should be an integer value");
                System.exit(1);
            }

            try {
                if ("null".equals(maxRandomString)) {
                    maxRandomString = "na";
                }
                maxRandom = Integer.parseInt(maxRandomString);
            } catch (NumberFormatException e) {
                System.err.println("Number of random interactions should be an integer value");
                System.exit(1);
            }

            if (!(minPathLength <= maxLength && maxPathLength <= maxLength && minPathLength <= maxPathLength && minPathLength > 0)) {
                System.err.println("Please check length of paths to report");
            }

            if ("null".equals(prefix)) {
                System.err.println("Prefix should be a string");
                System.exit(1);
            }

            String outputFile = of + File.separator + "masterPath_" + prefix.toUpperCase();
            String outputFolderRandom = outputFolder + File.separator + "Random" + File.separator;

            File randomFolder = new File(outputFolderRandom);
            randomFolder.mkdir();

            if ("integrated".equals(nwType)) {
                wr.find_shortest_paths_and_calculate_centrality(
                        nw,
                        hitList,
                        finalImplementers,
                        outputFile,
                        prefix,
                        maxLength,
                        minPathLength,
                        maxPathLength,
                        minCentrality,
                        outputFolderRandom,
                        prefix + "R",
                        maxRandom,
                        threads);
            } else {
                wr.find_shortest_paths_and_calculate_centrality_ppi(
                        nw,
                        hitList,
                        finalImplementers,
                        outputFile,
                        prefix,
                        maxLength,
                        minPathLength,
                        maxPathLength,
                        minCentrality,
                        outputFolderRandom,
                        prefix + "R",
                        maxRandom,
                        threads);
            }

        } else if ("network".equals(mode)) {
            File fi = new File(pathlist);

            if (!(fi.exists() && fi.isFile())) {
                System.err.println(pathlist);
                System.err.println("List of paths file does not exist. Mandatory for network mode");
                System.exit(1);
            }
            if ("null".equals(prefix)) {
                System.err.println("Prefix should be a string");
                System.exit(1);
            }
            String outputFile = of + File.separator + "masterPath_" + prefix.toUpperCase();
           
            
            putils.create_cyto_for_list_of_pathways(
                    outputFile + "_shortest_paths", 
                    pathlist, 
                    outputFolder, 
                    nw,
                    dbutils.hugo_by_id);
            
        }

        /*
           Network all0_hippie_h = new Network();
        Network all0_hippie_m0 = new Network();
        Network all0_hprd = new Network();
        Network all0_hippie_m = new Network();
        Network all0_hippie_meq = new Network();
        Network all0_hippie_a = new Network();
        Network all2 = new Network();
        Network all_invivo = new Network();
        
         //PPI network with high confidence interactions
         all0_hippie_h.integrateNetworks(nutils, putils, outils, dbutils);
         all0_hippie_h.saveNetworktofile(FoldersPaths.dbpath + "nodes_hippie_h", FoldersPaths.dbpath + "int_hippie_h");
         all0_hippie_h.loadNetworkfromfile(FoldersPaths.dbpath + "nodes_hippie_h", FoldersPaths.dbpath + "int_hippie_h");
         nutils.build_map_of_neighbors(all0_hippie_h);
       
         //PPI network with interactions for which confidence score is more 0
         all0_hippie_m0.integrateNetworks(nutils, putils, outils, dbutils);
         all0_hippie_m0.saveNetworktofile(FoldersPaths.dbpath + "nodes_hippie_m0", FoldersPaths.dbpath + "int_hippie_m0");
         all0_hippie_m0.loadNetworkfromfile(FoldersPaths.dbpath + "nodes_hippie_m0", FoldersPaths.dbpath + "int_hippie_m0");
         nutils.build_map_of_neighbors(all0_hippie_m0);
         
         //PPI network with all available interactions
         all0_hippie_a.integrateNetworks(nutils, putils, outils, dbutils);
         all0_hippie_a.saveNetworktofile(FoldersPaths.dbpath + "nodes_hippie_a", FoldersPaths.dbpath + "int_hippie_a");
         all0_hippie_a.loadNetworkfromfile(FoldersPaths.dbpath + "nodes_hippie_a", FoldersPaths.dbpath + "int_hippie_a");
         nutils.build_map_of_neighbors(all0_hippie_a);
      
         //PPI network with high + medium confidence interactions        
         all0_hippie_meq.integrateNetworks(nutils, putils, outils, dbutils);
         all0_hippie_meq.saveNetworktofile(FoldersPaths.dbpath + "nodes_hippie_meq", FoldersPaths.dbpath + "int_hippie_meq");
         all0_hippie_meq.loadNetworkfromfile(FoldersPaths.dbpath + "nodes_hippie_meq", FoldersPaths.dbpath + "int_hippie_meq");
         nutils.build_map_of_neighbors(all0_hippie_meq);
    
         //PPI network with only HPRD interactions        
         all0_hprd.integrateNetworks(nutils, putils, outils, dbutils);
         all0_hprd.saveNetworktofile(FoldersPaths.dbpath + "nodes_hprd", FoldersPaths.dbpath + "int_hprd");
         all0_hprd.loadNetworkfromfile(FoldersPaths.dbpath + "nodes_hprd", FoldersPaths.dbpath + "int_hprd");
         nutils.build_map_of_neighbors(all0_hprd);
         */
        //Integrated network        
        //all2.integrateNetworks(nutils, putils, outils, dbutils);
        //all2.saveNetworktofile(FoldersPaths.dbpath + "all_nodes_all_hippie_high_ogg1", FoldersPaths.dbpath + "all_int_all_hippie_high_ogg1");
        //all2.loadNetworkfromfile(FoldersPaths.dbpath + "all_nodes_all_hippie_high_ogg1", FoldersPaths.dbpath + "all_int_all_hippie_high_ogg1");
        //nutils.build_map_of_neighbors(all2);
        //In vivo ppi network        
        //all_invivo.integrateNetworks(nutils, putils, outils, dbutils);
        //all_invivo.saveNetworktofile(FoldersPaths.dbpath + "all_nodes_invivo", FoldersPaths.dbpath + "all_int_invivo");
        //all_invivo.loadNetworkfromfile(FoldersPaths.dbpath + "all_nodes_invivo", FoldersPaths.dbpath + "all_int_invivo");
        //nutils.build_map_of_neighbors(all_invivo);
    }

}
