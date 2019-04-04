package masterPATH;


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Network class describes an Network object
 * 
 * @author Natalia Rubanova
 */
public class Network {

    Map<String, Node> nodes;
    Map<String, Interaction> interactions;

    /**
     * Constructor
     */
    public Network() {
        this.nodes = new HashMap();
        this.interactions = new HashMap();
    }


    private void print_interactions_to_file(String file) throws IOException {
        //      String id,List<String> other_ids, Node int1, Node int2, String type, String sourcedb,List<String []> sourcedbentry, String quality, String dir
        BufferedWriter wr = new BufferedWriter(new FileWriter(file));
        for (String id : interactions.keySet()) {
            wr.write(id + "\t");
            for (String s : interactions.get(id).other_ids) {
                wr.write(s + ";");
            }

            //       wr.write(interactions.get(id).sourcedbentry.get(0)[0] + "\t" + interactions.get(i).sourcedbentry.get(0)[2] + "\t" + interactions.get(i).sourcedbentry.get(0)[5] + "\t" + "");
            wr.newLine();
        }
        wr.close();
    }

    /**
     * Load network from two files     *
     * @param nodef file with information about nodes
     * @param intf file with information about interactions
     */
    public void loadNetworkfromfile(String nodef, String intf) throws FileNotFoundException, IOException {
        //String id, String type, String id_type, String db_flag, List<String[]> ids  
        BufferedReader rn = new BufferedReader(new FileReader(nodef));
        BufferedReader ri = new BufferedReader(new FileReader(intf));
        String s, t;
        String[] ss, tt;
        Node n;
        Interaction in;
        List<String[]> ls;
        List<String> ls1, ls2;

        while ((s = rn.readLine()) != null) {
            ss = s.split("\t");
            ls = new ArrayList();
            try {
                tt = ss[4].split(";");
                for (String r : tt) {
                    ls.add(r.split("_"));
                }
                n = new Node(ss[0], ss[1], ss[2], ss[3], ls);
                this.nodes.put(ss[0], n);
            } catch (ArrayIndexOutOfBoundsException e) {
                //System.out.println(s);
                ls.add(new String[]{"_"});
                n = new Node(ss[0], ss[1], ss[2], ss[3], ls);
                this.nodes.put(ss[0], n);
            }

        }

        //String id, List<String> other_ids, Node int1, Node int2, String type, String sourcedb, 
        //List<String> sourcedbentry, String quality, String dir
        while ((s = ri.readLine()) != null) {
            ss = s.split("\t");
            ls2 = new ArrayList();
            ls1 = new ArrayList();
            tt = ss[1].split(";");
            ls1.addAll(Arrays.asList(tt));
            tt = ss[6].split(";");
            ls2.addAll(Arrays.asList(tt));
            in = new Interaction(ss[0], ls1, this.nodes.get(ss[2]), this.nodes.get(ss[3]), ss[4], ss[5], ls2, ss[7], ss[8]);
            this.interactions.put(ss[0], in);
        }
        rn.close();
        ri.close();
    }

    /**
     * Save all interactions and nodes in the network to the two output files
     * @param nodef node file
     * @param intf interaction file
     */
    public void saveNetworktofile(String nodef, String intf) throws IOException {
        BufferedWriter wn = new BufferedWriter(new FileWriter(nodef));
        BufferedWriter wi = new BufferedWriter(new FileWriter(intf));

        for (Node n : this.nodes.values()) {
            try {
                wn.write(n.toString());
            } catch (NullPointerException e) {
                System.out.println(n.id + n.db_flag);
            }
            wn.newLine();
        }

        for (Interaction in : this.interactions.values()) {
            wi.write(in.toString());
            wi.newLine();
        }

        wn.close();
        wi.close();
    }

    /**
     * Remove all interactions and nodes from network
     */
    public void removeAll() {
        for (String s : this.interactions.keySet()) {
            this.interactions.remove(s);
        }
        for (String s : this.nodes.keySet()) {
            this.nodes.remove(s);
        }
    }

    /**
     * Load network from databases
     *
     * @param nutils NetworkManager object
     * @param putils PathwayManager object
     * @param outils CentralityManager object
     * @param dbutils DBManager object
     */
    public void integrateNetworks(NetworkManager nutils, PathwayManager putils, CentralityManager outils, DBManager dbutils, int [] db) throws IOException {
        Network all=new Network();
        List<Network> networks = new ArrayList();
        
        if (db[0] == 1){
            dbutils.loadHPRD();
            networks.add(dbutils.hprd);
        }
        
        
        if (db[1] == 1) {
            dbutils.loadHIPPIE(false);
            networks.add(dbutils.hippie_high);
        }
        
        
        if (db[2] == 1) {
            dbutils.loadSignor();
            networks.add(dbutils.signor);
        }
        
        if (db[3] == 1) {
            dbutils.loadSignalink();
            networks.add(dbutils.signalink);
        }
        
        if (db[4] == 1) {
            dbutils.loadKEGG();
            networks.add(dbutils.kegg);
        }
        
        if(db[5] == 1){
            dbutils.loadTransmir();
            networks.add(dbutils.transmir);
        }
        
        if(db[6] == 1){
            dbutils.loadmirTarBase();
            networks.add(dbutils.mirtarbase);
        }
        
        if(db[7] == 1){
            dbutils.loadtFacts();
            networks.add(dbutils.tfacts);
        }
         
        
        
        
        
        //
        
         
       // networks.add(dbutils.hippie_medium_meq);
        // networks.add(dbutils.hippie_medium_m);
        //networks.add(dbutils.hippie_medium_m0);
       // networks.add(dbutils.hippie_invivo);
        
        
        //
        //
        //
        //
        //
        //
        
        all = nutils.merge_list_of_networks(networks);
        nutils.build_map_of_neighbors(all);
        nutils.add_missing_genes_and_products(all);
        //nutils.build_map_of_neighbors(all);
        this.interactions = all.interactions;
        this.nodes = all.nodes;
        all = null;

    }

}
