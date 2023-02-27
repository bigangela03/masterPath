/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package masterPATH;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author nrubanov
 */
public class RunnablePermuteLabelPPI {
    
    RandomManager rand = new RandomManager();
    Network nw;
    NetworkManager nutils; 
    PathwayManager putils; 
    CentralityManager outils; 
    DBManager dbutils; 
    String folder_for_random_paths; 
    int thread;
    String prefix;
    int number_of_permutations; 
    List<String> hit_list; 
    Map<String, String[]> hugo_by_id; 
    String file_fimplementers;
    int max_length_for_shortest_path;
    
    public RunnablePermuteLabelPPI( 
            Network nw, 
            NetworkManager nutils, 
            PathwayManager putils, 
            CentralityManager outils, 
            DBManager dbutils, 
            String folder_for_random_paths, 
            int thread,
            String prefix, 
            int number_of_permutations, 
            List<String> hit_list, 
            Map<String, String[]> hugo_by_id, 
            String file_fimplementers,
            int max_length_for_shortest_path) {
        
        this.nw = nw;
        this.nutils = nutils;
        this.putils = putils;
        this.outils = outils;
        this.dbutils = dbutils;
        this.folder_for_random_paths = folder_for_random_paths;
        this.thread = thread;
        this.prefix = prefix;
        this.number_of_permutations = number_of_permutations;
        this.hit_list = hit_list;
        this.hugo_by_id = hugo_by_id;
        this.file_fimplementers = file_fimplementers;
        this.max_length_for_shortest_path = max_length_for_shortest_path;
    }

    public void run()  {
        try {
            rand.permute_phenotype_label(
                    this.nw,
                    this.nutils,
                    this.putils,
                    this.outils,
                    this.dbutils,
                    this.folder_for_random_paths,
                    this.thread,
                    this.prefix ,
                    this.number_of_permutations,
                    this.hit_list,
                    this.dbutils.hugo_by_id,
                    this.file_fimplementers,
                    this.max_length_for_shortest_path
            );
        } catch (IOException ex) {
            Logger.getLogger(RunnablePermuteLabelPPI.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

}
