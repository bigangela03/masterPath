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
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

/**
 * NetworkManager class contains methods to deal with networks
 *
 * @author Natalia Rubanova
 */
public class NetworkManager {

    static Map<String, String> fpl = new HashMap();
    static Map<String, String> hg = new HashMap();

    /**
     * Merge list of networks
     *
     * @param n List of networks
     */
    public Network merge_list_of_networks(List<Network> n) {
        System.out.println("\n +++++++++++++merge list of networks++++++++++++ \n");
        Network tmp_add = new Network(), tmp_sub = new Network(), mrg = new Network();

        boolean match = false;
        int i = 0, j = 0, k = 0, l = 0;
        for (Network nw : n) {
            mrg.nodes.putAll(nw.nodes);
        }

        for (Network nw : n) {
            // tmp_add = new Network();
            // tmp_sub = new Network();
            for (Interaction int2 : nw.interactions.values()) {
                tmp_add.removeAll();
                tmp_sub.removeAll();
                for (Interaction int1 : mrg.interactions.values()) {
                    if (!int1.type.contains("PPI") || !int2.type.contains("PPI")) {
                        continue;
                    }
                    if ((int1.dir.equals("direct")) && (int2.dir.equals("direct"))) {
                        if ((int1.int1.id.equals(int2.int1.id) && int1.int2.id.equals(int2.int2.id))) {
                            match = true;
                            i++;
                            break;
                        }
                    } else {
                        if ((int1.int1.id.equals(int2.int1.id) && int1.int2.id.equals(int2.int2.id)) || (int1.int1.id.equals(int2.int2.id) && int1.int2.id.equals(int2.int1.id))) {
                            match = true;
                            i++;
                            if ((int1.dir.equals("indirect")) && (int2.dir.equals("direct"))) {
                                tmp_sub.interactions.put(int1.id, int1);
                                tmp_add.interactions.put(int2.id, int2);
                            }
                            if ((int1.dir.equals("direct")) && (int2.dir.equals("indirect"))) {
                                tmp_sub.interactions.put(int2.id, int2);
                                tmp_add.interactions.put(int1.id, int1);
                            }
                            break;
                        }
                    }
                }
                if (!match) {
                    tmp_add.interactions.put(int2.id, int2);
                }
                match = false;
                mrg.interactions.putAll(tmp_add.interactions);
                for (String int_id : tmp_sub.interactions.keySet()) {
                    mrg.interactions.remove(int_id);
                }
            }

        }

        for (Network nw : n) {
            System.out.println("network nodes " + nw.nodes.size() + " interactions " + nw.interactions.size());
        }
        System.out.println("Overlap interactions " + i);
        return mrg;
    }

    /**
     * Create links to all neighbors for every node in the network
     *
     * @param nw Network
     */
    public void build_map_of_neighbors(Network nw) {
        System.out.println("+++++++++++buildNeighbours++++++++++");
        Map<String, String> it_type = new HashMap();
        String n1;
        String n2;
        Map<String, Interaction> tmp_map = new HashMap();
        List<String> duplicate_int = new ArrayList();
        int sum_degree = 0, sum_degree_up = 0, sum_degree_down = 0, sum_degree_rev = 0, indirect_count = 0, direct_count = 0;
        for (String s : nw.nodes.keySet()) {
            sum_degree = sum_degree + nw.nodes.get(s).downnbrs.size()
                    + nw.nodes.get(s).upnbrs.size()
                    + nw.nodes.get(s).revnbrs.size();
            nw.nodes.get(s).downnbrs = new HashMap();
            nw.nodes.get(s).upnbrs = new HashMap();
            nw.nodes.get(s).revnbrs = new HashMap();
        }
        System.out.println("sum degree1 " + sum_degree);
        System.out.println("interactions1 " + nw.interactions.size());

        for (String it : nw.interactions.keySet()) {
            if (nw.interactions.get(it).int1.type.toLowerCase().equals("stimulus")
                    || nw.interactions.get(it).int1.type.toLowerCase().equals("chemical")
                    || nw.interactions.get(it).int1.type.toLowerCase().equals("smallmolecule")
                    || nw.interactions.get(it).int1.type.toLowerCase().equals("phenotype")
                    || nw.interactions.get(it).int1.type.toLowerCase().equals("proteinfamily")
                    || nw.interactions.get(it).int2.type.toLowerCase().equals("stimulus")
                    || nw.interactions.get(it).int2.type.toLowerCase().equals("chemical")
                    || nw.interactions.get(it).int2.type.toLowerCase().equals("smallmolecule")
                    || nw.interactions.get(it).int2.type.toLowerCase().equals("phenotype")
                    || nw.interactions.get(it).int2.type.toLowerCase().equals("proteinfamily")) {
                //duplicate_int.add(it);
                continue;
            }

            if (!nw.interactions.get(it).int1.type.toLowerCase().equals("stimulus")
                    && !nw.interactions.get(it).int1.type.toLowerCase().equals("chemical")
                    && !nw.interactions.get(it).int1.type.toLowerCase().equals("smallmolecule")
                    && !nw.interactions.get(it).int1.type.toLowerCase().equals("phenotype")
                    && !nw.interactions.get(it).int1.type.toLowerCase().equals("proteinfamily")
                    && !nw.interactions.get(it).int2.type.toLowerCase().equals("stimulus")
                    && !nw.interactions.get(it).int2.type.toLowerCase().equals("chemical")
                    && !nw.interactions.get(it).int2.type.toLowerCase().equals("smallmolecule")
                    && !nw.interactions.get(it).int2.type.toLowerCase().equals("phenotype")
                    && !nw.interactions.get(it).int2.type.toLowerCase().equals("proteinfamily")
                    && !nw.interactions.get(it).int2.type.toLowerCase().equals("protein")
                    && !nw.interactions.get(it).int2.type.toLowerCase().equals("gene")
                    && !nw.interactions.get(it).int2.type.toLowerCase().equals("compound")
                    && !nw.interactions.get(it).int2.type.toLowerCase().equals("glycan")
                    && !nw.interactions.get(it).int2.type.toLowerCase().equals("mirna")
                    && !nw.interactions.get(it).int2.type.toLowerCase().equals("complex")
                    && !nw.interactions.get(it).int1.type.toLowerCase().equals("protein")
                    && !nw.interactions.get(it).int1.type.toLowerCase().equals("gene")
                    && !nw.interactions.get(it).int1.type.toLowerCase().equals("compound")
                    && !nw.interactions.get(it).int1.type.toLowerCase().equals("glycan")
                    && !nw.interactions.get(it).int1.type.toLowerCase().equals("mirna")
                    && !nw.interactions.get(it).int1.type.toLowerCase().equals("complex")) {
                System.out.println(nw.interactions.get(it).int1.type + " " + nw.interactions.get(it).int2.type);
                break;
            }
            try {
                n1 = nw.interactions.get(it).int1.id;
                n2 = nw.interactions.get(it).int2.id;

            } catch (NullPointerException e) {
                System.out.print("BN error " + it + "\t");
                System.out.println(nw.interactions.get(it).sourcedbentry.get(0));
                continue;
            }
            if (n1.equals(n2)) {
                if (!duplicate_int.contains(it)) {
                    //duplicate_int.add(it);
                }
                continue;
            }

            if (nw.interactions.get(it).dir.toLowerCase().equals("indirect") || nw.interactions.get(it).type.toLowerCase().equals("kegg_reaction_reversible")) {

                if (nw.nodes.get(n1).revnbrs.containsKey(n2)) {
                    if (!duplicate_int.contains(it)) {
                        //   duplicate_int.add(it);
                    }
                    continue;
                }

                if (nw.nodes.get(n2).revnbrs.containsKey(n1)) {
                    if (!duplicate_int.contains(it)) {
                        //   duplicate_int.add(it);
                    }
                    continue;
                }
                indirect_count++;
                if (nw.nodes.get(n1).revnbrs.containsKey(n2) || nw.nodes.get(n2).revnbrs.containsKey(n1)) {
                    System.out.println("ERROR2 ");
                    break;
                }

                nw.nodes.get(n1).revnbrs.put(n2, nw.interactions.get(it));
                nw.nodes.get(n2).revnbrs.put(n1, nw.interactions.get(it));

            } else {

                if (nw.nodes.get(n1).downnbrs.containsKey(n2)) {
                    //  System.out.println("AAAA11 " + nw.nodes.get(n1).downnbrs.get(nw.interactions.get(it).int2.id).type);
                    //  System.out.println("AAAA12 " + nw.interactions.get(it).type);
                    if (!duplicate_int.contains(it)) {
                        //   duplicate_int.add(it);
                    }
                    continue;
                }
                if (nw.nodes.get(n2).upnbrs.containsKey(n1)) {
                    //  System.out.println("AAAA21 " + nw.nodes.get(n1).upnbrs.get(nw.interactions.get(it).int1.id).id);
                    //  System.out.println("AAAA22 " + it);
                    if (!duplicate_int.contains(it)) {
                        //   duplicate_int.add(it);
                    }
                    continue;
                }

                if (nw.nodes.get(n1).downnbrs.containsKey(n2) && !nw.nodes.get(n2).upnbrs.containsKey(n1)) {
                    System.out.println("ERROR ");
                    break;
                }
                direct_count++;
                if (nw.nodes.get(n1).downnbrs.containsKey(n2) || nw.nodes.get(n2).upnbrs.containsKey(n1)) {
                    System.out.println("ERROR3 ");
                    break;
                }
                nw.nodes.get(n1).downnbrs.put(n2, nw.interactions.get(it));
                nw.nodes.get(n2).upnbrs.put(n1, nw.interactions.get(it));

            }
            /*if (nw.interactions.get(it).type.equals("kegg_reaction_reversible")) {
             //nw.nodes.get(n1).revnbrs.put(nw.interactions.get(it).int2, nw.interactions.get(it));
             nw.nodes.get(n1).revnbrs.put(nw.interactions.get(it).int2.id, nw.interactions.get(it));
             //nw.nodes.get(n2).revnbrs.put(nw.interactions.get(it).int1, nw.interactions.get(it));
             nw.nodes.get(n2).revnbrs.put(nw.interactions.get(it).int1.id, nw.interactions.get(it));
             } else {
             //nw.nodes.get(n1).downnbrs.put(nw.interactions.get(it).int2, nw.interactions.get(it));
             nw.nodes.get(n1).downnbrs.put(nw.interactions.get(it).int2.id, nw.interactions.get(it));
             //nw.nodes.get(n2).upnbrs.put(nw.interactions.get(it).int1, nw.interactions.get(it));
             nw.nodes.get(n2).upnbrs.put(nw.interactions.get(it).int1.id, nw.interactions.get(it));
             }*/
        }
        for (String s : duplicate_int) {
            nw.interactions.remove(s);
        }
        for (String s : it_type.keySet()) {
            // System.out.println(s);
        }
        sum_degree = 0;
        for (String s : nw.nodes.keySet()) {
            sum_degree = sum_degree + nw.nodes.get(s).downnbrs.size()
                    + nw.nodes.get(s).upnbrs.size()
                    + nw.nodes.get(s).revnbrs.size();
            sum_degree_down = sum_degree_down + nw.nodes.get(s).downnbrs.size();
            sum_degree_up = sum_degree_up + nw.nodes.get(s).upnbrs.size();
            sum_degree_rev = sum_degree_rev + nw.nodes.get(s).revnbrs.size();
        }
        System.out.println("indirect count " + indirect_count + " direct_count " + direct_count);
        System.out.println("sum degree2 " + sum_degree);
        System.out.println("sum degree down " + sum_degree_down);
        System.out.println("sum degree up " + sum_degree_up);
        System.out.println("sum degree rev " + sum_degree_rev);
        System.out.println("interactions2 " + nw.interactions.size());
    }

    /**
     * Load hit genes and final implementers files
     *
     * @param hgf Hit genes file
     * @param fpf Final players file
     * @param hugo Link to the HGNC nomenclature
     */
    public void load_hitlist_and_finalimpl(String hgf, String fpf, Map<String, String[]> hugo) throws FileNotFoundException, IOException {
        System.out.println(" +++++++++++++loadHitGenes_and_FinalPlayers++++++++++++ ");
        System.out.println("Loading hit list and fplayers files");
        char[] alphabet = "123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ".toCharArray();
        char[] digit = "123456789".toCharArray();
        //loadHitGenes_and_FinalPlayers(hgf, fpf, hugo);
        fpl = new HashMap<String, String>();
        hg = new HashMap<String, String>();
        String s;
        String[] ss;
        int i = 0;
        int j = 0;
        BufferedReader rd = new BufferedReader(new FileReader(fpf));
        while ((s = rd.readLine()) != null) {
            ss = s.split("\t");
            fpl.put(ss[1], ss[0]);
        }
        rd.close();
        rd = new BufferedReader(new FileReader(hgf));
        boolean flag = false;
        String sep = "";
        while ((s = rd.readLine()) != null) {
            i++;
            if (hugo.containsKey(s)) {
                flag = false;
                if (s.startsWith("MIR") && (hugo.get(s)[2].startsWith("microRNA"))) {
                    if (s.contains("LET")) {
                        hg.put("hsa-let-" + s.substring(6).toLowerCase(), s);
                        flag = true;
                        //   System.out.println(s + "\t" + "hsa-let-" + s.substring(6));
                    } else {
                        flag = true;
                        hg.put("hsa-mir-" + s.substring(3).toLowerCase(), s);
                        //  System.out.println(s + "\t" + "hsa-mir-" + s.substring(3));
                    }
                } else {
                    flag = true;
                    hg.put("p" + hugo.get(s)[0], "p" + s);// hit gene - gene
                }
            } else if (s.startsWith("MIR")) {
                flag = false;
                for (char c : alphabet) {
                    if (s.contains("LET")) {
                        if (hugo.containsKey(s + c) && !(Character.isDigit(c) && Character.isDigit(s.charAt(s.length() - 1)))) {
                            if (Character.isDigit(c)) {
                                hg.put("hsa-let-" + s.substring(6).toLowerCase() + "-" + c, s + c);
                            } else {
                                hg.put("hsa-let-" + s.substring(6).toLowerCase() + "" + c, s + c);
                            }
                            flag = true;
                        }
                        if (hugo.containsKey(s + "-" + c)) {
                            hg.put("hsa-let-" + s.substring(6).toLowerCase() + "-" + c, s);
                            flag = true;
                        }
                    } else {
                        if (hugo.containsKey(s + c) && !(Character.isDigit(c) && Character.isDigit(s.charAt(s.length() - 1)))) {
                            if (Character.isDigit(c)) {
                                hg.put("hsa-mir-" + s.substring(3).toLowerCase() + "-" + c, s + c);
                            } else {
                                hg.put("hsa-mir-" + s.substring(3).toLowerCase() + c, s + c);
                            }
                            flag = true;
                        }
                        if (hugo.containsKey(s + "-" + c)) {
                            hg.put("hsa-mir-" + s.substring(3).toLowerCase(), s);
                            flag = true;
                        }
                    }
                }
            }
            if (!flag) {
                System.out.println("Not in HUGO! " + s);
            } else {
                j++;
            }

        }
        /*if (s.equals("MIR206")) {
         hg.put("hsa-mir-206", s);
         } */

        rd.close();

        System.out.println(j
                + "/" + i + " finished");
    }

    private int return_node_type(String current_node) {
        int node_type;
        //System.out.println(current_node);
        if (current_node.toLowerCase().equals("protein") || current_node.toLowerCase().equals("complex") || current_node.toLowerCase().equals("proteinfamily")) {
            node_type = 0;
        } else if (current_node.toLowerCase().equals("gene")) {
            node_type = 1;
        } else if (current_node.toLowerCase().equals("mirna")) {
            node_type = 2;
        } else {
            node_type = 3;
        }
        return node_type;
    }

    /**
     * Find length-bound pathways between two lists of nodes
     *
     * @param network Network
     * @param n1 List of Hit genes nodes
     * @param n2 List of Final implementers nodes
     * @param max Maximum length
     * @param prefix Prefix for pathway id
     * @param out Output file
     */
    public List<PathUnit> find_pathway_BF_algorithm(Network network, List<String> n1, List<String> n2, int max, String prefix, FileWriter out) throws IOException {
        PathUnit curn, tmp_curn;//Buffered
        List<String> path, seeds;
        LinkedList<PathUnit> found = new LinkedList(), pre_found = new LinkedList();

        int current_hg = 0, id = 0;
        Map<String, Integer> visited_depth = new HashMap();
        String tmp_s, path_string = "";
        List<PathUnit> tmp_paths, paths = new ArrayList();
        List<String> vn, tmp_vn;
        Map<String, int[]> visited_nodes, tmp_visited_nodes;
        int[] tmp_path_type;
        int[] only_ppi = new int[]{1, 0, 0, 0};
        int path_type;
        boolean cont;
        int continue_length;
        int jj;
        List<String> fruitfull_seeds = new ArrayList();
        for (String beg : n1) {
            //pre_found=new LinkedList();
            visited_nodes = new HashMap();
            current_hg++;
            System.out.println("Current hit gene " + beg + "(" + current_hg + "/" + n1.size() + ")");
            continue_length = 0;
////////////////////////////////////////////////////////////////////////////////
            if (network.nodes.containsKey(beg)) {
                for (String n : network.nodes.get(beg).downnbrs.keySet()) {
                    path = new ArrayList();
                    path.add(network.nodes.get(beg).downnbrs.get(n).id);
                    // path_string=network.nodes.get(beg).downnbrs.get(n).id;
                    seeds = new ArrayList();
                    seeds.add(n);
                    tmp_curn = new PathUnit(path, seeds, n, network.nodes.get(beg).downnbrs.get(n).id);
                    paths.add(tmp_curn);
                    if (n2.contains(n)) {
                        pre_found.add(tmp_curn);
                    }
                    path_type = return_node_type(network.nodes.get(n).type);
                    if ((!visited_nodes.containsKey(n))) {
                        tmp_path_type = new int[]{0, 0, 0, 0};
                        tmp_path_type[path_type] = 1;
                        visited_nodes.put(n, tmp_path_type);
                    }

                }
////////////////////////////////////////////////////////////////////////////////
                for (String n : network.nodes.get(beg).revnbrs.keySet()) {
                    path = new ArrayList();
                    path.add(network.nodes.get(beg).revnbrs.get(n).id);
                    // path_string=network.nodes.get(beg).revnbrs.get(n).id;
                    seeds = new ArrayList();
                    seeds.add(n);
                    tmp_curn = new PathUnit(path, seeds, n, network.nodes.get(beg).revnbrs.get(n).id);
                    paths.add(tmp_curn);
                    if (n2.contains(n)) {
                        pre_found.add(tmp_curn);
                    }
                    path_type = return_node_type(network.nodes.get(n).type);
                    if ((!visited_nodes.containsKey(n))) {
                        tmp_path_type = new int[]{0, 0, 0, 0};
                        tmp_path_type[path_type] = 1;
                        visited_nodes.put(n, tmp_path_type);
                    }
                }
            }
////////////////////////////////////////////////////////////////////////////////
            for (int i = 2; i <= 5; i++) {
                int p_s = paths.size();
                System.out.println("length " + i + "\t" + "paths size " + paths.size());
                tmp_paths = new ArrayList();
                tmp_visited_nodes = new HashMap();
                for (int j = 0; j < p_s; j++) {

                    curn = paths.get(j);
////////////////////////////////////////////////////////////////////////////////           
                    for (String n : network.nodes.get(curn.current_seed).revnbrs.keySet()) {
                        // path_string="";
                        if (!network.nodes.get(curn.current_seed).revnbrs.get(n).id.equals(curn.path.get(curn.path.size() - 1))) {
                            path = new ArrayList();
                            path.addAll(curn.path);
                            path.add(network.nodes.get(curn.current_seed).revnbrs.get(n).id);
                            path_string = curn.path_string + "\t" + network.nodes.get(curn.current_seed).revnbrs.get(n).id;
                            seeds = new ArrayList();
                            seeds.addAll(curn.seeds);
                            seeds.add(n);
                            tmp_curn = new PathUnit(path, seeds, n, path_string);
                            path_type = return_node_type(network.nodes.get(n).type);
                            if (n2.contains(n)) {
                                pre_found.add(tmp_curn);
                            }
                            if (curn.seeds.contains(n)) {
                                continue;
                            }
                            cont = true;
                            if (visited_nodes.containsKey(n)) {

                                if (!Arrays.equals(only_ppi, visited_nodes.get(curn.current_seed))) {
                                    cont = false;
                                }
                                if (cont) {
                                    continue;
                                }
                            }
                            if (visited_depth.containsKey(n) && (visited_depth.get(n) + i > max)) {
                                continue_length++;
                                continue;
                            }
                            if (i < 5) {
                                tmp_paths.add(tmp_curn);
                                if (visited_nodes.containsKey(n) && Arrays.equals(only_ppi, visited_nodes.get(curn.current_seed))) {
                                    tmp_path_type = visited_nodes.get(n).clone();
                                } else {
                                    tmp_path_type = visited_nodes.get(curn.current_seed).clone();
                                }
                                tmp_path_type[path_type] = 1;
                                if (!Arrays.equals(only_ppi, tmp_path_type) || !tmp_visited_nodes.containsKey(n)) {
                                    tmp_visited_nodes.put(n, tmp_path_type);
                                }
                            }
                        }
                    }
////////////////////////////////////////////////////////////////////////////////
                    for (String n : network.nodes.get(curn.current_seed).downnbrs.keySet()) {
                        //path_string="";
                        path = new ArrayList();
                        path.addAll(curn.path);
                        path.add(network.nodes.get(curn.current_seed).downnbrs.get(n).id);
                        path_string = curn.path_string + "\t" + network.nodes.get(curn.current_seed).downnbrs.get(n).id;
                        seeds = new ArrayList();
                        seeds.addAll(curn.seeds);
                        seeds.add(n);
                        tmp_curn = new PathUnit(path, seeds, n, path_string);
                        path_type = return_node_type(network.nodes.get(n).type);
                        if (n2.contains(n)) {
                            pre_found.add(tmp_curn);
                        }
                        if (curn.seeds.contains(n)) {
                            continue;
                        }
                        cont = true;
                        if (visited_nodes.containsKey(n)) {//!Arrays.equals(visited_nodes.get(n), visited_nodes.get(curn.current_seed))  || || path_type >= 1 
                            if (!Arrays.equals(only_ppi, visited_nodes.get(curn.current_seed))) {
                                cont = false;
                            }
                            if (cont) {
                                continue;
                            }
                        }
                        if (visited_depth.containsKey(n) && (visited_depth.get(n) + i > max)) {
                            continue_length++;
                            continue;
                        }
                        if (i < 5) {
                            tmp_paths.add(tmp_curn);
                            if (visited_nodes.containsKey(n) && Arrays.equals(only_ppi, visited_nodes.get(curn.current_seed))) {
                                tmp_path_type = visited_nodes.get(n).clone();
                            } else {
                                tmp_path_type = visited_nodes.get(curn.current_seed).clone();
                            }
                            tmp_path_type[path_type] = 1;
                            //!tmp_visited_nodes.containsKey(n)
                            if (!Arrays.equals(only_ppi, tmp_path_type) || !tmp_visited_nodes.containsKey(n)) {
                                tmp_visited_nodes.put(n, tmp_path_type);

                            }
                        }
                    }
                    if (paths.size() > 50000) {
                        if ((j % 50000) == 0) {
                            System.out.print(j + "\t");
                        }
                    }
                }
                //System.out.println("continue_length " + continue_length);
                if (i < 5) {
                    paths = new ArrayList();
                    paths.addAll(tmp_paths);
                }
                visited_nodes.putAll(tmp_visited_nodes);
            }
////////////////////////////////////////////////////////////////////////////////
            System.out.println("Fruitfull seeds");
            for (PathUnit pu : pre_found) {
                for (String seed : pu.seeds) {
                    if (!fruitfull_seeds.contains(seed)) {
                        fruitfull_seeds.add(seed);
                    }
                    visited_depth.remove(seed);
                }
            }
            System.out.println("Visited depth");
            if (paths.size() <= 200000) {
                for (PathUnit pu : paths) {
                    jj = 0;
                    for (String seed : pu.seeds) {
                        if (!fruitfull_seeds.contains(seed) && (!visited_depth.containsKey(seed) || (visited_depth.containsKey(seed) && (visited_depth.get(seed) < pu.seeds.size() - jj - 1)))) {
                            visited_depth.put(seed, pu.seeds.size() - jj - 1);
                        }
                        jj++;
                    }
                }
            }
            System.out.println("Writing" + pre_found.size());

            //          for (int i = 0; i < pre_found.size(); i++) {
            for (PathUnit pu : pre_found) {
                if ((id % 10000) == 0) {
                    System.out.println(id);
                }
                id++;
                //tmp_s = prefix + id + "\t" + beg;
                //  for (int j = 0; j < pre_found.get(i).path.size(); j++) {
                //      tmp_s = tmp_s + "\t" + pre_found.get(i).path.get(j);
                //  }
                //  for (String t : pre_found.get(i).path) {
                //      tmp_s = tmp_s + "\t" + t;
                //  }
                tmp_s = prefix + id + "\t" + beg + "\t" + pu.path_string + "\t" + pu.current_seed + "\n";
                out.write(tmp_s);
                //  out.newLine();
            }
            System.out.println("New vars");
            found = new LinkedList();
            pre_found = new LinkedList();
            paths = new ArrayList();
        }
        path = null;
        System.gc();
        return pre_found;
    }

    /**
     * Find length-bound pathways between two lists of nodes in a PPI network
     *
     * @param network
     * @param n1 List of Hit genes nodes
     * @param n2 List of Final implementers nodes
     * @param max Maximum length
     * @param prefix Prefix for interaction id
     * @param out output file
     */
    public List<PathUnit> find_pathway_BF_algorithm_ppi(Network network, List<String> n1, List<String> n2, int max, String prefix, BufferedWriter out) throws IOException {
        PathUnit curn, tmp_curn;
        List<String> path, seeds;
        LinkedList<PathUnit> found = new LinkedList(), pre_found = new LinkedList();

        int current_hg = 0, id = 0;
        Map<String, Integer> visited_depth = new HashMap();
        String tmp_s, path_string = "";
        List<PathUnit> tmp_paths, paths = new ArrayList();
        List<String> vn, tmp_vn;
        Map<String, int[]> visited_nodes, tmp_visited_nodes;
        int[] tmp_path_type;
        int[] only_ppi = new int[]{1, 0, 0, 0};
        int path_type;
        boolean cont;
        int continue_length;
        int jj;
        List<String> fruitfull_seeds = new ArrayList();
        for (String beg : n1) {
            n2.remove(beg);
            if (n2.isEmpty()) {
                break;
            }
            visited_nodes = new HashMap();
            current_hg++;
            System.out.println("Current hit gene " + beg + "(" + current_hg + "/" + n1.size() + ")");
            continue_length = 0;
////////////////////////////////////////////////////////////////////////////////
            if (network.nodes.containsKey(beg)) {
                for (String n : network.nodes.get(beg).downnbrs.keySet()) {
                    path = new ArrayList();
                    path.add(network.nodes.get(beg).downnbrs.get(n).id);
                    seeds = new ArrayList();
                    seeds.add(n);
                    tmp_curn = new PathUnit(path, seeds, n, network.nodes.get(beg).downnbrs.get(n).id);
                    paths.add(tmp_curn);
                    if (n2.contains(n)) {
                        pre_found.add(tmp_curn);
                    }
                    path_type = return_node_type(network.nodes.get(n).type);
                    if ((!visited_nodes.containsKey(n))) {
                        tmp_path_type = new int[]{0, 0, 0, 0};
                        tmp_path_type[path_type] = 1;
                        visited_nodes.put(n, tmp_path_type);
                    }

                }
////////////////////////////////////////////////////////////////////////////////
                for (String n : network.nodes.get(beg).revnbrs.keySet()) {
                    path = new ArrayList();
                    path.add(network.nodes.get(beg).revnbrs.get(n).id);
                    seeds = new ArrayList();
                    seeds.add(n);
                    tmp_curn = new PathUnit(path, seeds, n, network.nodes.get(beg).revnbrs.get(n).id);
                    paths.add(tmp_curn);
                    if (n2.contains(n)) {
                        pre_found.add(tmp_curn);
                    }
                    path_type = return_node_type(network.nodes.get(n).type);
                    if ((!visited_nodes.containsKey(n))) {
                        tmp_path_type = new int[]{0, 0, 0, 0};
                        tmp_path_type[path_type] = 1;
                        visited_nodes.put(n, tmp_path_type);
                    }
                }
            }
////////////////////////////////////////////////////////////////////////////////
            for (int i = 2; i <= 5; i++) {
                int p_s = paths.size();
                System.out.println("length " + i + "\t" + "paths size " + paths.size());
                tmp_paths = new ArrayList();
                tmp_visited_nodes = new HashMap();
                for (int j = 0; j < p_s; j++) {
                    curn = paths.get(j);
////////////////////////////////////////////////////////////////////////////////           
                    for (String n : network.nodes.get(curn.current_seed).revnbrs.keySet()) {
                        if (!network.nodes.get(curn.current_seed).revnbrs.get(n).id.equals(curn.path.get(curn.path.size() - 1))) {
                            path = new ArrayList();
                            path.addAll(curn.path);
                            path.add(network.nodes.get(curn.current_seed).revnbrs.get(n).id);
                            path_string = curn.path_string + "\t" + network.nodes.get(curn.current_seed).revnbrs.get(n).id;
                            seeds = new ArrayList();
                            seeds.addAll(curn.seeds);
                            seeds.add(n);
                            tmp_curn = new PathUnit(path, seeds, n, path_string);
                            path_type = return_node_type(network.nodes.get(n).type);
                            if (n2.contains(n)) {
                                pre_found.add(tmp_curn);
                            }
                            if (curn.seeds.contains(n)) {
                                continue;
                            }
                            cont = true;
                            if (visited_nodes.containsKey(n)) {//!Arrays.equals(visited_nodes.get(n), visited_nodes.get(curn.current_seed))  || || path_type >= 1 
                                if (!Arrays.equals(only_ppi, visited_nodes.get(curn.current_seed))) {
                                    cont = false;
                                }
                                if (cont) {
                                    continue;
                                }
                            }
                            if (visited_depth.containsKey(n) && (visited_depth.get(n) + i > max)) {
                                continue_length++;
                                continue;
                            }
                            if (i < 5) {
                                tmp_paths.add(tmp_curn);
                                if (visited_nodes.containsKey(n) && Arrays.equals(only_ppi, visited_nodes.get(curn.current_seed))) {
                                    tmp_path_type = visited_nodes.get(n).clone();
                                } else {
                                    tmp_path_type = visited_nodes.get(curn.current_seed).clone();
                                }
                                tmp_path_type[path_type] = 1;
                                //||  !tmp_visited_nodes.containsKey(n)
                                if (!Arrays.equals(only_ppi, tmp_path_type) || !tmp_visited_nodes.containsKey(n)) {
                                    tmp_visited_nodes.put(n, tmp_path_type);
                                }
                            }
                        }
                    }
////////////////////////////////////////////////////////////////////////////////
                    for (String n : network.nodes.get(curn.current_seed).downnbrs.keySet()) {
                        path = new ArrayList();
                        path.addAll(curn.path);
                        path.add(network.nodes.get(curn.current_seed).downnbrs.get(n).id);
                        path_string = curn.path_string + "\t" + network.nodes.get(curn.current_seed).downnbrs.get(n).id;
                        seeds = new ArrayList();
                        seeds.addAll(curn.seeds);
                        seeds.add(n);
                        tmp_curn = new PathUnit(path, seeds, n, path_string);
                        path_type = return_node_type(network.nodes.get(n).type);
                        if (n2.contains(n)) {
                            pre_found.add(tmp_curn);
                        }
                        if (curn.seeds.contains(n)) {
                            continue;
                        }
                        cont = true;
                        if (visited_nodes.containsKey(n)) {//!Arrays.equals(visited_nodes.get(n), visited_nodes.get(curn.current_seed))  || || path_type >= 1 
                            if (!Arrays.equals(only_ppi, visited_nodes.get(curn.current_seed))) {
                                cont = false;
                            }
                            if (cont) {
                                continue;
                            }
                        }
                        if (visited_depth.containsKey(n) && (visited_depth.get(n) + i > max)) {
                            continue_length++;
                            continue;
                        }
                        if (i < 5) {
                            tmp_paths.add(tmp_curn);
                            if (visited_nodes.containsKey(n) && Arrays.equals(only_ppi, visited_nodes.get(curn.current_seed))) {
                                tmp_path_type = visited_nodes.get(n).clone();
                            } else {
                                tmp_path_type = visited_nodes.get(curn.current_seed).clone();
                            }
                            tmp_path_type[path_type] = 1;
                            //!tmp_visited_nodes.containsKey(n)
                            if (!Arrays.equals(only_ppi, tmp_path_type) || !tmp_visited_nodes.containsKey(n)) {
                                tmp_visited_nodes.put(n, tmp_path_type);

                            }
                        }
                    }
                    if (paths.size() > 50000) {
                        if ((j % 50000) == 0) {
                            System.out.print(j + "\t");
                        }
                    }

                }
                //System.out.println("continue_length " + continue_length);
                if (i < 5) {
                    paths = new ArrayList();
                    paths.addAll(tmp_paths);
                }
                visited_nodes.putAll(tmp_visited_nodes);
            }
////////////////////////////////////////////////////////////////////////////////
            System.out.println("Fruitfull seeds");
            for (PathUnit pu : pre_found) {
                for (String seed : pu.seeds) {
                    if (!fruitfull_seeds.contains(seed)) {
                        fruitfull_seeds.add(seed);
                    }
                    visited_depth.remove(seed);
                }
            }
            System.out.println("Visited depth");
            if (paths.size() <= 200000) {
                for (PathUnit pu : paths) {
                    jj = 0;
                    for (String seed : pu.seeds) {
                        if (!fruitfull_seeds.contains(seed) && (!visited_depth.containsKey(seed) || (visited_depth.containsKey(seed) && (visited_depth.get(seed) < pu.seeds.size() - jj - 1)))) {
                            visited_depth.put(seed, pu.seeds.size() - jj - 1);
                        }
                        jj++;
                    }
                }
            }
            System.out.println("Writing" + pre_found.size());

            //          for (int i = 0; i < pre_found.size(); i++) {
            for (PathUnit pu : pre_found) {
                if ((id % 10000) == 0) {
                    System.out.println(id);
                }
                id++;
                //tmp_s = prefix + id + "\t" + beg;
                //  for (int j = 0; j < pre_found.get(i).path.size(); j++) {
                //      tmp_s = tmp_s + "\t" + pre_found.get(i).path.get(j);
                //  }
                //  for (String t : pre_found.get(i).path) {
                //      tmp_s = tmp_s + "\t" + t;
                //  }
                tmp_s = prefix + id + "\t" + beg + "\t" + pu.path_string + "\t" + pu.current_seed + "\n";
                out.write(tmp_s);
                //  out.newLine();
            }
            System.out.println("New vars");
            found = new LinkedList();
            pre_found = new LinkedList();
            paths = new ArrayList();
        }
        path = null;
        System.gc();
        return pre_found;
    }

    /**
     * Wrapper for finding pathways
     *
     * @param nw Network
     * @param hugo Link to HGNC nomenclature
     * @param max Maximum length
     * @param outf Output file
     * @param prefix Prefix for interaction id
     */
    public void find_pathway_for_list_BF_algorithm(Network nw, Map<String, String[]> hugo, int max, String outf, String prefix) throws IOException {
        System.out.println("++++++++++getPathforListBF+++++++++");
        FileWriter out = new FileWriter(outf);
        List<PathUnit> paths;
        int id = 0;
        int current_hg = 0;
        // for (String s : hg.keySet()) {
        //current_hg++;
        //System.out.println("Current hit gene " + s + "(" + current_hg + "/" + hg.keySet().size() + ")");
        paths = find_pathway_BF_algorithm(nw, new ArrayList(hg.keySet()), new ArrayList(fpl.keySet()), max, prefix, out);
        /*for (int i = 0; i < paths.size(); i++) {
         id++;
         out.write(prefix + id + "\t");
         out.write("s");
         for (int j = 0; j < paths.get(i).path.size(); j++) {
         out.write("\t" + paths.get(i).path.get(j));
         }
         out.write("\t" + paths.get(i).current_seed);
         out.newLine();
         }
         out.flush();
         // } */
        out.close();
    }

    /**
     * Wrapper for finding pathways in a PPI network
     *
     * @param nw Network
     * @param hugo Link to HGNC nomenclature
     * @param max Maximum length
     * @param outf Output file
     * @param prefix Prefix for interaction id
     */
    public void find_pathway_for_list_BF_algorithm_ppi(Network nw, Map<String, String[]> hugo, int max, String outf, String prefix) throws IOException {
        System.out.println("++++++++++getPathforListBF+++++++++");
        BufferedWriter out = new BufferedWriter(new FileWriter(outf));
        List<PathUnit> paths;
        int id = 0;
        int current_hg = 0;
        paths = find_pathway_BF_algorithm_ppi(nw, new ArrayList(hg.keySet()), new ArrayList(fpl.keySet()), max, prefix, out);
        out.close();
    }

    /**
     * Add transcription information
     *
     * @param nw Network
     */
    public void add_missing_genes_and_products(Network nw) {
        System.out.println("+++++++++++add_genes_and_products++++++++++");
        System.out.println("Before statistics: nodes " + nw.nodes.size() + " interactions " + nw.interactions.size());
        String id, type, id_type, tmp_id;
        Node n2;
        Interaction int1;
        int i = 0;
        List<String> other_ids;

        List<String> sourcedbentry;
        Map<String, Node> new_nodes = new HashMap();
        Map<String, Interaction> tmp_map;
        Map<String, Interaction> new_int = new HashMap();
        boolean found;
        Node n1;
        Node[] nn = new Node[2];
        List<String> lookedat = new ArrayList();
        for (Interaction tr : nw.interactions.values()) {
            nn[0] = tr.int1;
            nn[1] = tr.int2;
            for (Node n : nn) {
                found = false;
                id = n.id;
                type = n.type;
                id_type = n.id_type;
                if (type.toLowerCase().equals("protein") && !lookedat.contains(id)) {
                    //&& id_type.toLowerCase().contains("hugo")
                    tmp_id = id.substring(1);
                    if (!nw.nodes.containsKey(tmp_id)) {
                        i++;
                        n2 = new Node(tmp_id, "gene", "hugo", "10000008", n.ids);
                        other_ids = new ArrayList();
                        other_ids.add("-");
                        sourcedbentry = new ArrayList();
                        sourcedbentry.add("-");
                        int1 = new Interaction("man" + i, other_ids, n2, n, "expression", "man", sourcedbentry, "1", "direct");
                        n2.downnbrs.put(n.id, int1);
                        nw.nodes.get(n.id).upnbrs.put(n2.id, int1);
                        new_int.put("man" + i, int1);
                        new_nodes.put(tmp_id, n2);
                        lookedat.add(id);
                        lookedat.add(tmp_id);
                        found = true;
                        continue;
                    }
                    if (nw.nodes.containsKey(tmp_id)) {
                        for (String ss : n.upnbrs.keySet()) {
                            if (ss.equals(tmp_id)) {
                                found = true;
                                lookedat.add(id);
                                lookedat.add(tmp_id);
                                continue;
                            }
                        }
                        if (!found) {
                            i++;
                            other_ids = new ArrayList();
                            other_ids.add("-");
                            sourcedbentry = new ArrayList();
                            int1 = new Interaction("man" + i, other_ids, nw.nodes.get(tmp_id), n, "expression", "man", sourcedbentry, "1", "direct");
                            nw.nodes.get(tmp_id).downnbrs.put(n.id, int1);
                            nw.nodes.get(n.id).upnbrs.put(tmp_id, int1);
                            new_int.put("man" + i, int1);
                            lookedat.add(id);
                            lookedat.add(tmp_id);
                            continue;
                        }
                    }
                }
                if (type.toLowerCase().equals("gene") && !lookedat.contains(id)) {
                    //&& id_type.toLowerCase().contains("hugo")
                    tmp_id = "p" + id;
                    if (!nw.nodes.containsKey(tmp_id)) {
                        i++;
                        n2 = new Node(tmp_id, "protein", "hugo", "10000008", n.ids);
                        other_ids = new ArrayList();
                        other_ids.add("-");
                        sourcedbentry = new ArrayList();
                        int1 = new Interaction("man" + i, other_ids, n, n2, "expression", "man", sourcedbentry, "1", "direct");
                        n2.upnbrs.put(n.id, int1);
                        nw.nodes.get(n.id).downnbrs.put(n2.id, int1);
                        new_nodes.put(tmp_id, n2);
                        new_int.put("man" + i, int1);
                        lookedat.add(id);
                        lookedat.add(tmp_id);
                        found = true;
                        continue;
                    }
                    if (nw.nodes.containsKey(tmp_id)) {
                        for (String ss : n.downnbrs.keySet()) {
                            if (ss.equals(tmp_id)) {
                                found = true;
                                lookedat.add(id);
                                lookedat.add(tmp_id);
                                continue;
                            }
                        }
                        if (!found) {
                            i++;
                            other_ids = new ArrayList();
                            other_ids.add("-");
                            sourcedbentry = new ArrayList();
                            int1 = new Interaction("man" + i, other_ids, n, nw.nodes.get(tmp_id), "expression", "man", sourcedbentry, "1", "direct");
                            nw.nodes.get(tmp_id).upnbrs.put(n.id, int1);
                            nw.nodes.get(n.id).downnbrs.put(tmp_id, int1);
                            new_int.put("man" + i, int1);
                            lookedat.add(id);
                            lookedat.add(tmp_id);
                            continue;
                        }
                    }
                }
            }
        }

        nw.nodes.putAll(new_nodes);
        nw.interactions.putAll(new_int);
        System.out.println("After statistics: nodes " + nw.nodes.size() + " interactions " + nw.interactions.size());
        System.out.println("Added nodes " + new_nodes.size() + " int  " + new_int.size());
        System.out.println("i " + i);
    }

    private String db_name_by_number(int i) {
        switch (i) {
            case 2:
                return "HPRD";
            case 3:
                return "TP";
            case 4:
                return "KEGG";
            case 5:
                return "TRANSMIR";
            case 6:
                return "MIRTARBASE";
            case 7:
                return "TFACTS";
            default: {
                System.out.println("error");
                return "error";

            }
        }
    }

    class Depth_List_of_PathUnit {

        int max_depth;
        List<PathUnit> punits;

        public Depth_List_of_PathUnit() {
            this.max_depth = 0;
            this.punits = new ArrayList();
        }

        public Depth_List_of_PathUnit(int max_depth, List<PathUnit> punits) {
            this.max_depth = max_depth;
            this.punits = punits;
        }
    }

    class PathUnit {

        List<String> path;
        List<String> seeds; // list of seeds
        String current_seed;
        String path_string;

        public PathUnit() {
            this.path = new ArrayList();
            this.seeds = new ArrayList();
            this.current_seed = "";
        }

        public PathUnit(List<String> path, List<String> seeds, String seed, String path_string) {
            this.path = path;
            this.current_seed = seed;
            this.seeds = seeds;
            this.path_string = path_string;
        }

        public String toString() {
            String s = this.current_seed;
            for (int i = 0; i < this.path.size(); i++) {
                s = s + "\t" + this.path.get(i);
            }
            return s;
        }
    }

    class PathUnitwFP {

        List<String> path;
        String seed;
        String fp;

        public PathUnitwFP() {
            this.path = new ArrayList();
            this.seed = "";
            this.fp = "";
        }

        public PathUnitwFP(List<String> path, String seed, String fp) {
            this.path = path;
            this.seed = seed;
            this.fp = fp;
        }

        public String toString() {
            String s = this.seed;
            for (int i = 0; i < this.path.size(); i++) {
                s = s + "\t" + this.path.get(i);
            }
            s = s + "\t" + this.fp;
            return s;
        }
    }
}

/*
 private void copy_files(String f1, String f2) throws IOException {
 String s;
 BufferedReader source = new BufferedReader(new FileReader(f1));
 BufferedWriter target = new BufferedWriter(new FileWriter(f2));
 while ((s = source.readLine()) != null) {
 target.write(s);
 target.newLine();
 }
 source.close();
 target.close();
 }


 for (Node n : nw.nodes.get("pHGNC:7611").downrevnbrs.keySet()) {
 /// System.out.println(n.id);
 }


 m1 = merged.interactions.get(m_id).int1;
 m2 = merged.interactions.get(m_id).int2;
 try {
 boolean c = db1.id.equals(db2.id);
 // boolean d= m1.id.equals(m2.id);
 // System.out.println();
 } catch (NullPointerException e) {
 //  System.out.println("m_id=  " + m_id);
 System.out.println("id=  " + id);
 //   break;
 }
 try {
 boolean d = m1.id.equals(m2.id);
 // System.out.println();
 } catch (NullPointerException e) {
 System.out.println("m_id=  " + m_id);
 // break;
 }






 try {n2 = nw.interactions.get(it).int2.id;}
 catch(NullPointerException e){
 System.out.println("wwwwwwwwwwwww " +n1+"  "+ it);
 break;
 }
 */
/*
 public List<List<String>> getPath(Network nw, String n1, String n2, int max) {
 System.out.println("++++++++++getPath+++++++++");

 List<String> path = new ArrayList();
 List<String> tmp = new ArrayList();
 List<List<String>> paths = new ArrayList();// interactions
 List<List<String>> tmps = new ArrayList(); // interactions
 List<List<String>> found = new ArrayList(); // interactions
 List<List<String>> seeds= new ArrayList(); // 
 List<String> s_tmp =new ArrayList();

 for (Node n : nw.nodes.get(n1).downrevnbrs.keySet()) {
 tmp = new ArrayList();
 tmp.add(nw.nodes.get(n1).downrevnbrs.get(n).id);
 paths.add(tmp);
 s_tmp=new ArrayList();
 s_tmp.add(n1);
 seeds.add(s_tmp);
 }
 //System.out.println(" -1- " + paths.size());

 Node lastn, prev;
 String prevprev;
 for (int i = 1; i <= max; i++) {
 //System.out.println("  length " + i + "  paths.size() " + paths.size());
 tmps = new ArrayList();
 tmps.addAll(paths);
 //System.out.println("            tmps.size() " + tmps.size());
 paths = new ArrayList();
 for (int j = 0; j < tmps.size(); j++) {
                
 lastn = nw.interactions.get(tmps.get(j).get(i - 1)).int2;
 prev = nw.interactions.get(tmps.get(j).get(i - 1)).int1;
                
 if (i > 1) {
 prevprev = nw.interactions.get(tmps.get(j).get(i - 2)).int2.id;
 } else {
 prevprev = "!";
 }
 // System.out.println(" int1 " + nw.interactions.get(tmps.get(j).get(i-1 )).int2.id);
 if (lastn.id.equals(n2)) {
 tmp = new ArrayList();
 tmp.addAll(tmps.get(j));
 found.add(tmp);
 } else {
 //  System.out.println("            lastn.downrevnbrs.size    " + nw.nodes.get(lastn.id).downrevnbrs.size());
 for (Node n : nw.nodes.get(lastn.id).downrevnbrs.keySet()) {
 if ((!prev.id.equals(nw.nodes.get(lastn.id).downrevnbrs.get(n).int2.id)) && (!prevprev.equals(nw.nodes.get(lastn.id).downrevnbrs.get(n).int2.id)) && (!prevprev.equals(nw.nodes.get(lastn.id).downrevnbrs.get(n).int1.id))) {
 tmp = new ArrayList();
 tmp.addAll(tmps.get(j));
 // System.out.println("            tmp.size    " + tmp.size());
 tmp.add(nw.nodes.get(lastn.id).downrevnbrs.get(n).id);
 paths.add(tmp);
 }
 }
 }
 }
 }

 return found;
 }

 */
/**
 *
 * public void getLongPathBFd(Network nw, String n1, String n2, int max, String
 * tempfolder, BufferedWriter found_out) throws IOException { //
 * System.out.println("++++++++++getLongPath+++++++++"); Node lastn, prev;
 * String prevprev; List<String> tmp = new ArrayList(); List<PathUnit> paths =
 * new ArrayList();// interactions List<PathUnit> tmps = new ArrayList(); //
 * interactions List<PathUnit> found = new ArrayList(); // interactions PathUnit
 * pu; String pathsf = tempfolder + "paths.txt"; String tmpsf = tempfolder +
 * "tmps.txt"; String tmp_s; BufferedWriter paths_out; BufferedWriter tmps_out;
 * BufferedReader paths_in; BufferedReader tmps_in; BufferedReader found_in;
 *
 * paths_out = new BufferedWriter(new FileWriter(pathsf)); /*for (Node n :
 * nw.nodes.get(n1).downrevnbrs.keySet()) { //tmp = new ArrayList();
 * //tmp.add(nw.nodes.get(n1).downrevnbrs.get(n).id); //paths.add(new
 * PathUnit(tmp, n1)); // pu = new PathUnit(tmp, n1); //
 * paths_out.write(pu.toString()); paths_out.write(n1 + "\t" +
 * nw.nodes.get(n1).downrevnbrs.get(n).id); paths_out.newLine(); }
 */
/*    for (Node n : nw.nodes.get(n1).downnbrs.keySet()) {
 //tmp = new ArrayList();
 //tmp.add(nw.nodes.get(n1).downrevnbrs.get(n).id);
 //paths.add(new PathUnit(tmp, n1));           
 // pu = new PathUnit(tmp, n1);
 // paths_out.write(pu.toString());
 paths_out.write(n1 + "\t" + nw.nodes.get(n1).downnbrs.get(n).id);
 paths_out.newLine();
 }
 for (Node n : nw.nodes.get(n1).revnbrs.keySet()) {
 //tmp = new ArrayList();
 //tmp.add(nw.nodes.get(n1).downrevnbrs.get(n).id);
 //paths.add(new PathUnit(tmp, n1));           
 // pu = new PathUnit(tmp, n1);
 // paths_out.write(pu.toString());
 paths_out.write(n1 + "\t" + nw.nodes.get(n1).revnbrs.get(n).id);
 paths_out.newLine();
 }

 paths_out.close();
 for (int i = 1; i <= max; i++) {
 //tmps = new ArrayList();
 //tmps.addAll(paths);
 copy_files(pathsf, tmpsf);
 // paths = new ArrayList();
 paths_out = new BufferedWriter(new FileWriter(pathsf));
 tmps_in = new BufferedReader(new FileReader(tmpsf));
 String s;
 String[] ss;
 int j = 0;
 //for (int j = 0; j < tmps.size(); j++) {
 while ((s = tmps_in.readLine()) != null) {
 ss = s.split("\t");
 // System.out.println(s);
 //lastn = nw.interactions.get(tmps.get(j).path.get(i - 1)).int2;                
 lastn = nw.interactions.get(ss[i - 1 + 1]).int2;

 //if (tmps.get(j).seed.equals(lastn.id)) {
 if (ss[0].equals(lastn.id)) {
 prev = lastn;
 //lastn = nw.interactions.get(tmps.get(j).path.get(i - 1)).int1;
 lastn = nw.interactions.get(ss[i - 1 + 1]).int1;
 } else {
 //prev = nw.interactions.get(tmps.get(j).path.get(i - 1)).int1;
 prev = nw.interactions.get(ss[i - 1 + 1]).int1;
 }
 if (i > 1) {
 //prevprev = nw.interactions.get(tmps.get(j).path.get(i - 2)).int2.id;
 prevprev = nw.interactions.get(ss[i - 2 + 1]).int2.id;
 } else {
 prevprev = "!";
 }
 if (lastn.id.equals(n2)) {
 //tmp = new ArrayList();
 tmp_s = ss[1];
 //tmp.addAll(tmps.get(j).path);
 for (int k = 2; k < ss.length; k++) {
 tmp_s = tmp_s + "\t" + ss[k];
 }
 //found.add(new PathUnit(tmp, lastn.id));
 found_out.write(lastn.id + "\t" + tmp_s);
 found_out.newLine();
 } else {
 /*for (Node n : nw.nodes.get(lastn.id).downrevnbrs.keySet()) {
 if ((!prev.id.equals(nw.nodes.get(lastn.id).downrevnbrs.get(n).int2.id)) && (!prevprev.equals(nw.nodes.get(lastn.id).downrevnbrs.get(n).int2.id)) && (!prevprev.equals(nw.nodes.get(lastn.id).downrevnbrs.get(n).int1.id))) {
 //tmp = new ArrayList();
 tmp_s = ss[1];
 //tmp.addAll(tmps.get(j).path);
 for (int k = 2; k < ss.length; k++) {
 tmp_s = tmp_s + "\t" + ss[k];
 }
 //tmp.add(nw.nodes.get(lastn.id).downrevnbrs.get(n).id);
 tmp_s = tmp_s + "\t" + nw.nodes.get(lastn.id).downrevnbrs.get(n).id;
 //paths.add(new PathUnit(tmp, lastn.id));
 paths_out.write(lastn.id + "\t" + tmp_s);
 paths_out.newLine();
 }
 }*/
/*     for (Node n : nw.nodes.get(lastn.id).revnbrs.keySet()) {
 if ((!prev.id.equals(nw.nodes.get(lastn.id).revnbrs.get(n).int2.id)) && (!prevprev.equals(nw.nodes.get(lastn.id).revnbrs.get(n).int2.id)) && (!prevprev.equals(nw.nodes.get(lastn.id).revnbrs.get(n).int1.id))) {
 tmp_s = ss[1];
 for (int k = 2; k < ss.length; k++) {
 tmp_s = tmp_s + "\t" + ss[k];
 }
 tmp_s = tmp_s + "\t" + nw.nodes.get(lastn.id).revnbrs.get(n).id;
 paths_out.write(lastn.id + "\t" + tmp_s);
 paths_out.newLine();
 }
 }
 for (Node n : nw.nodes.get(lastn.id).downnbrs.keySet()) {
 tmp_s = ss[1];
 for (int k = 2; k < ss.length; k++) {
 tmp_s = tmp_s + "\t" + ss[k];
 }
 tmp_s = tmp_s + "\t" + nw.nodes.get(lastn.id).downnbrs.get(n).id;
 paths_out.write(lastn.id + "\t" + tmp_s);
 paths_out.newLine();
 }
 }
 }
 paths_out.close();
 }
 // return found;
 }

 public void getLongPathsforListBFd(Network nw, String hgf, String fpf, Map<String, String[]> hugo, int max, String outfolder) throws IOException {
 System.out.println("++++++++++getLongPathforLIst+++++++++");
 System.out.println("Loading hit list and fplayers files");
 load_hitlist_and_finalimpl(hgf, fpf, hugo);
 System.out.println("finished");
 String foundf = outfolder + "found.txt";
 BufferedWriter found_out = new BufferedWriter(new FileWriter(foundf));
 for (String s : fpl.keySet()) {
 for (String t : hg.keySet()) {
 System.out.println(fpl.get(s) + "   " + hg.get(t));
 getLongPathBFd(nw, s, t, max, outfolder, found_out);
 }
 }
 found_out.close();
 }


 public List<PathUnit> getPath(Network nw, String n1, String n2, int max) {
 System.out.println("++++++++++getPath+++++++++");
 List<String> tmp = new ArrayList();
 List<PathUnit> paths = new ArrayList();// interactions
 List<PathUnit> tmps = new ArrayList(); // interactions
 List<PathUnit> found = new ArrayList(); // interactions

 /*for (Node n : nw.nodes.get(n1).downrevnbrs.keySet()) {
 tmp = new ArrayList();
 tmp.add(nw.nodes.get(n1).downrevnbrs.get(n).id);
 paths.add(new PathUnit(tmp, n1));
 }*/
/*      for (Node n : nw.nodes.get(n1).downnbrs.keySet()) {
 tmp = new ArrayList();
 tmp.add(nw.nodes.get(n1).downnbrs.get(n).id);
 paths.add(new PathUnit(tmp, n1));
 }
 for (Node n : nw.nodes.get(n1).revnbrs.keySet()) {
 tmp = new ArrayList();
 tmp.add(nw.nodes.get(n1).revnbrs.get(n).id);
 paths.add(new PathUnit(tmp, n1));
 }
 // System.out.println(" -1- " + paths.size() + " " + n1 + " " + nw.nodes.get(n1).downrevnbrs.size()+ " " + n2 +" "+ nw.nodes.get(n2).downrevnbrs.size());

 Node lastn, prev;
 String prevprev;
 for (int i = 1; i <= max; i++) {
 tmps = new ArrayList();
 tmps.addAll(paths);
 paths = new ArrayList();
 for (int j = 0; j < tmps.size(); j++) {
 lastn = nw.interactions.get(tmps.get(j).path.get(i - 1)).int2;
 if (tmps.get(j).seed.equals(lastn.id)) {
 prev = lastn;
 lastn = nw.interactions.get(tmps.get(j).path.get(i - 1)).int1;
 } else {
 prev = nw.interactions.get(tmps.get(j).path.get(i - 1)).int1;
 }
 if (i > 1) {
 prevprev = nw.interactions.get(tmps.get(j).path.get(i - 2)).int2.id;
 } else {
 prevprev = "!";
 }
 if (lastn.id.equals(n2)) {
 tmp = new ArrayList();
 tmp.addAll(tmps.get(j).path);
 found.add(new PathUnit(tmp, lastn.id));
 } else {
 /* for (Node n : nw.nodes.get(lastn.id).downrevnbrs.keySet()) {
 if ((!prev.id.equals(nw.nodes.get(lastn.id).downrevnbrs.get(n).int2.id)) && (!prevprev.equals(nw.nodes.get(lastn.id).downrevnbrs.get(n).int2.id)) && (!prevprev.equals(nw.nodes.get(lastn.id).downrevnbrs.get(n).int1.id))) {
 tmp = new ArrayList();
 tmp.addAll(tmps.get(j).path);
 // System.out.println("            tmp.size    " + tmp.size());
 tmp.add(nw.nodes.get(lastn.id).downrevnbrs.get(n).id);
 paths.add(new PathUnit(tmp, lastn.id));
 }
 }*/
/*               for (Node n : nw.nodes.get(lastn.id).revnbrs.keySet()) {
 if ((!prev.id.equals(nw.nodes.get(lastn.id).revnbrs.get(n).int2.id)) && (!prevprev.equals(nw.nodes.get(lastn.id).revnbrs.get(n).int2.id)) && (!prevprev.equals(nw.nodes.get(lastn.id).revnbrs.get(n).int1.id))) {
 tmp = new ArrayList();
 tmp.addAll(tmps.get(j).path);
 // System.out.println("            tmp.size    " + tmp.size());
 tmp.add(nw.nodes.get(lastn.id).revnbrs.get(n).id);
 paths.add(new PathUnit(tmp, lastn.id));
 }
 }
 for (Node n : nw.nodes.get(lastn.id).downnbrs.keySet()) {
 tmp = new ArrayList();
 tmp.addAll(tmps.get(j).path);
 // System.out.println("            tmp.size    " + tmp.size());
 tmp.add(nw.nodes.get(lastn.id).downnbrs.get(n).id);
 paths.add(new PathUnit(tmp, lastn.id));
 }
 }
 }
 }
 return found;
 }

 public void getPathsforList(Network nw, String hgf, String fpf, Map<String, String[]> hugo, int max, String outf) throws IOException {
 System.out.println("++++++++++getPathforLIst+++++++++");
 BufferedWriter out = new BufferedWriter(new FileWriter(outf));
 System.out.println("Loading hit list and fplayers files");
 load_hitlist_and_finalimpl(hgf, fpf, hugo);
 System.out.println("finished");
 List<PathUnit> paths = new ArrayList();
 for (String s : fpl.keySet()) {
 for (String t : hg.keySet()) {
 System.out.println(fpl.get(s) + "   " + hg.get(t));
 paths = getPath(nw, s, t, max);
 System.out.println(paths.size());
 for (int i = 0; i < paths.size(); i++) {
 out.write(s + " begin");
 for (int j = 0; j < paths.get(i).path.size(); j++) {
 out.write("\t" + nw.interactions.get(paths.get(i).path.get(j)).int1.id + "-" + nw.interactions.get(paths.get(i).path.get(j)).int2.id);
 //out.write(nw.interactions.get(paths.get(i).get(j)).sourcedbentry.get(0)[1] + "\t" + nw.interactions.get(paths.get(i).get(j)).sourcedbentry.get(0)[1] + "-" + nw.interactions.get(paths.get(i).get(j)).sourcedbentry.get(0)[4] + "\t-------");
 }
 out.write("\t + end " + t);
 out.newLine();
 }
 }
 }
 out.close();
 }



 public List<PathUnit> getLongPathDF(Network nw, String n1, String n2, int max) {
 // System.out.println("++++++++++getLongPathDF+++++++++");
 List<String> path;
 Stack<PathUnit> stack = new Stack();
 Stack<PathUnit> found = new Stack();
 PathUnit curn;
 //System.out.println("dwn " + nw.nodes.get(n1).downrevnbrs.size());
 //System.out.println("up " + nw.nodes.get(n1).upnbrs.size());
 /*for (Node n : nw.nodes.get(n1).downrevnbrs.keySet()) {
 path = new ArrayList();
 path.add(nw.nodes.get(n1).downrevnbrs.get(n).id);
 stack.push(new PathUnit(path, n.id));
 }*/
/*       if (nw.nodes.containsKey(n1)) {
 for (Node n : nw.nodes.get(n1).downnbrs.keySet()) {
 path = new ArrayList();
 path.add(nw.nodes.get(n1).downnbrs.get(n).id);
 stack.push(new PathUnit(path, n.id));
 }
 for (Node n : nw.nodes.get(n1).revnbrs.keySet()) {
 path = new ArrayList();
 path.add(nw.nodes.get(n1).revnbrs.get(n).id);
 stack.push(new PathUnit(path, n.id));
 }
 }
 while (!stack.empty()) {
 // System.out.println("Stack size " + stack.size());
 curn = stack.pop();
 if (curn.seed.equals(n2)) {
 found.push(curn);
 continue;
 }
 if (curn.path.size() >= max) {
 continue;
 }
 /*for (Node n : nw.nodes.get(curn.seed).downrevnbrs.keySet()) {
 if (!nw.nodes.get(curn.seed).downrevnbrs.get(n).id.equals(curn.path.get(curn.path.size() - 1))) {
 path = new ArrayList();
 path.addAll(curn.path);
 path.add(nw.nodes.get(curn.seed).downrevnbrs.get(n).id);
 stack.push(new PathUnit(path, n.id));
 }
 } */
/*           for (Node n : nw.nodes.get(curn.seed).revnbrs.keySet()) {
 if (!nw.nodes.get(curn.seed).revnbrs.get(n).id.equals(curn.path.get(curn.path.size() - 1))) {
 path = new ArrayList();
 path.addAll(curn.path);
 path.add(nw.nodes.get(curn.seed).revnbrs.get(n).id);
 stack.push(new PathUnit(path, n.id));
 }
 }
 for (Node n : nw.nodes.get(curn.seed).downnbrs.keySet()) {
 path = new ArrayList();
 path.addAll(curn.path);
 path.add(nw.nodes.get(curn.seed).downnbrs.get(n).id);
 stack.push(new PathUnit(path, n.id));
 }
 }
 stack=null;
 path=null;
 System.gc();
 return found;
 }



 public void getLongPathsforListDF(Network nw, String hgf, String fpf, Map<String, String[]> hugo, int max, String outf) throws IOException {
 System.out.println("++++++++++getPathforListDF+++++++++");
 BufferedWriter out = new BufferedWriter(new FileWriter(outf));
 System.out.println("Loading hit list and fplayers files");
 load_hitlist_and_finalimpl(hgf, fpf, hugo);
 System.out.println("finished");
 List<PathUnit> paths = new ArrayList();
 for (String t : fpl.keySet()) {
 for (String s : hg.keySet()) {
 // System.out.println(hg.get(s) + "   " + fpl.get(t));
 paths = getLongPathDF(nw, s, t, max);
 //  System.out.println(paths.size());
 for (int i = 0; i < paths.size(); i++) {
 //out.write(s + " begin");
 out.write(s);
 for (int j = 0; j < paths.get(i).path.size(); j++) {
 out.write("\t" + paths.get(i).path.get(j));
 }
 //out.write("\t + end " + t);
 out.write("\t" + t);
 out.newLine();
 }
 }
 }
 out.close();
 }

 public void load_hitlist_and_finalimpl(String hgf, String fpf, Map<String, String[]> hugo) throws FileNotFoundException, IOException {
 String s;
 String[] ss = new String[3];
 BufferedReader rd = new BufferedReader(new FileReader(fpf));

 while ((s = rd.readLine()) != null) {
 ss = s.split("\t");
 fpl.put(ss[1], ss[0]);
 }
 rd.close();

 //BufferedWriter wr = new BufferedWriter(new FileWriter(hgf+"hugo"));
 rd = new BufferedReader(new FileReader(hgf));
 //String hugo_id;
 while ((s = rd.readLine()) != null) {
 //hugo_id="-";
 if (hugo.containsKey(s)) {
 //hg.put("p" + hugo.get(s)[0], s); //hit gene - protein
 hg.put(hugo.get(s)[0], s);// hit gene - gene
 //hugo_id="p" + hugo.get(s)[0];
 }
 //  wr.write(s + "\t"+ hugo_id);
 //  wr.newLine();
 }
 rd.close();
 // wr.close();

 }



 public void add_missing_genes_and_products(Network nw) {
 System.out.println("+++++++++++add_missing_genes_and_products++++++++++");
 System.out.println("Before statistics: nodes " + nw.nodes.size() + " interactions " + nw.interactions.size());
 String id, type, id_type, gp_id;
 Node n2;
 int i = 0;
 List<String> other_ids;
 List<String[]> sourcedbentry;
 Map<String, Node> new_nodes = new HashMap();
 Map<String, Interaction> new_int = new HashMap();
 for (Node n : nw.nodes.values()) {
 id = n.id;
 type = n.type;
 id_type = n.id_type;
 if (type.toLowerCase().equals("protein")) {
 if (id_type.toLowerCase().contains("hugo")) {
 gp_id = id.substring(1);
 } else {
 gp_id = id;
 }
 if (!nw.nodes.containsKey(gp_id)) {
 i++;
 n2 = new Node(gp_id, "gene", "hugo", "10000008", n.ids);
 new_nodes.put(gp_id, n2);
 //nw.nodes.put(gp_id, n2);
 other_ids = new ArrayList();
 sourcedbentry = new ArrayList();
 new_int.put("man" + i, new Interaction("man" + i, other_ids, n2, n, "expression", "man", sourcedbentry, "manexpr", "-"));
 //nw.interactions.put("man" + i, new Interaction("man" + i, other_ids, n2, n, "expression", "man", sourcedbentry, "manexpr", "-"));
 } else {
 if (!nw.nodes.get(id).upnbrs.containsKey(nw.nodes.get(gp_id))) {
 i++;
 n2 = nw.nodes.get(gp_id);
 other_ids = new ArrayList();
 sourcedbentry = new ArrayList();
 //    new_int.put("man" + i, new Interaction("man" + i, other_ids, n2, n, "expression", "man", sourcedbentry, "manexpr", "-"));
 }
 }
 } else if (n.type.toLowerCase().equals("gene")) {
 if (id_type.toLowerCase().contains("hugo")) {
 gp_id = "p" + id;
 if (!nw.nodes.containsKey(gp_id)) {
 i++;
 n2 = new Node(gp_id, "protein", "hugo", "10000008", n.ids);
 new_nodes.put(gp_id, n2);
 //nw.nodes.put(gp_id, n2);
 other_ids = new ArrayList();
 sourcedbentry = new ArrayList();
 new_int.put("man" + i, new Interaction("man" + i, other_ids, n, n2, "expression", "man", sourcedbentry, "manexpr", "-"));
 //nw.interactions.put("man" + i, new Interaction("man" + i, other_ids, n, n2, "expression", "man", sourcedbentry, "manexpr", "-"));
 } else {
 if (!nw.nodes.get(id).upnbrs.containsKey(nw.nodes.get(gp_id))) {
 i++;
 n2 = nw.nodes.get(gp_id);
 other_ids = new ArrayList();
 sourcedbentry = new ArrayList();
 // new_int.put("man" + i, new Interaction("man" + i, other_ids, n, n2, "expression", "man", sourcedbentry, "manexpr", "-"));
 }
 }
 }
 }
 }
 nw.nodes.putAll(new_nodes);
 nw.interactions.putAll(new_int);
 System.out.println("After statistics: nodes " + nw.nodes.size() + " interactions " + nw.interactions.size());
 System.out.println("Added nodes " + new_nodes.size() + " int  " + new_int.size());
 System.out.println("i " + i);
 }

 */
/*
 for (Node n : nw.nodes.values()) {
 found = false;
 id = n.id;
 type = n.type;
 id_type = n.id_type;
 if (type.toLowerCase().equals("protein") && id_type.toLowerCase().contains("hugo")) {
 tmp_id = id.substring(1);
 if (!nw.nodes.containsKey(tmp_id)) {
 i++;
 n2 = new Node(tmp_id, "gene", "hugo", "10000008", n.ids);
 other_ids = new ArrayList();
 sourcedbentry = new ArrayList();
 int1 = new Interaction("man" + i, other_ids, n2, n, "expression", "man", sourcedbentry, "manexpr", "-");
 n2.downnbrs.put(n.id, int1);
 nw.nodes.get(n.id).upnbrs.put(n2.id, int1);
 new_int.put("man" + i, int1);
 new_nodes.put(tmp_id, n2);
 //new_int.put("man" + i, new Interaction("man" + i, other_ids, n2, n, "expression", "man", sourcedbentry, "manexpr", "-"));
 }
 if (nw.nodes.containsKey(tmp_id)) {
 tmp_map = new HashMap();
 for (String ss : n.upnbrs.keySet()) {
 if (ss.equals(tmp_id)) {
 found = true;
 break;
 }
 if (!found) {
 i++;
 other_ids = new ArrayList();
 sourcedbentry = new ArrayList();
 int1 = new Interaction("man" + i, other_ids, nw.nodes.get(tmp_id), n, "expression", "man", sourcedbentry, "manexpr", "-");
 // n.upnbrs.put(tmp_id, int1);
 // n2 = nw.nodes.get(tmp_id);
 // n2.downnbrs.put(n.id, int1);
 nw.nodes.get(tmp_id).downnbrs.put(n.id, int1);
 tmp_map.put(tmp_id, int1);
 //nw.nodes.get(n.id).upnbrs.put(tmp_id, int1);
 new_int.put("man" + i, int1);
 }
 }
 nw.nodes.get(n.id).upnbrs.putAll(tmp_map);
 tmp_map = null;
 }
 }
 if (type.toLowerCase().equals("gene") && id_type.toLowerCase().contains("hugo")) {
 tmp_id = "p" + id;
 if (!nw.nodes.containsKey(tmp_id)) {
 i++;
 n2 = new Node(tmp_id, "protein", "hugo", "10000008", n.ids);
 other_ids = new ArrayList();
 sourcedbentry = new ArrayList();
 int1 = new Interaction("man" + i, other_ids, n, n2, "expression", "man", sourcedbentry, "manexpr", "-");
 n2.upnbrs.put(n.id, int1);
 nw.nodes.get(n.id).downnbrs.put(n2.id, int1);
 new_nodes.put(tmp_id, n2);
 new_int.put("man" + i, int1);
 }
 if (nw.nodes.containsKey(tmp_id)) {
 tmp_map = new HashMap();
 for (String ss : n.downnbrs.keySet()) {
 if (ss.equals(tmp_id)) {
 found = true;
 break;
 }
 if (!found) {
 i++;
 other_ids = new ArrayList();
 sourcedbentry = new ArrayList();
 int1 = new Interaction("man" + i, other_ids, n, nw.nodes.get(tmp_id), "expression", "man", sourcedbentry, "manexpr", "-");
 nw.nodes.get(tmp_id).upnbrs.put(n.id, int1);
 //nw.nodes.get(n.id).downnbrs.put(tmp_id, int1);
 tmp_map.put(tmp_id, int1);
 new_int.put("man" + i, int1);
 }
 }
 nw.nodes.get(n.id).downnbrs.putAll(tmp_map);
 tmp_map = null;
 }
 }

 /* if (type.toLowerCase().equals("protein")) {
 if (id_type.toLowerCase().contains("hugo")) {
 gp_id = id.substring(1);
 } else {
 gp_id = id;
 }
 if (!nw.nodes.containsKey(gp_id)) {
 i++;
 n2 = new Node(gp_id, "gene", "hugo", "10000008", n.ids);
 new_nodes.put(gp_id, n2);
 //nw.nodes.put(gp_id, n2);
 other_ids = new ArrayList();
 sourcedbentry = new ArrayList();
 new_int.put("man" + i, new Interaction("man" + i, other_ids, n2, n, "expression", "man", sourcedbentry, "manexpr", "-"));
 //nw.interactions.put("man" + i, new Interaction("man" + i, other_ids, n2, n, "expression", "man", sourcedbentry, "manexpr", "-"));
 } else {
 if (!nw.nodes.get(id).upnbrs.containsKey(nw.nodes.get(gp_id))) {
 i++;
 n2 = nw.nodes.get(gp_id);
 other_ids = new ArrayList();
 sourcedbentry = new ArrayList();
 //    new_int.put("man" + i, new Interaction("man" + i, other_ids, n2, n, "expression", "man", sourcedbentry, "manexpr", "-"));
 }
 }
 } else if (n.type.toLowerCase().equals("gene")) {
 if (id_type.toLowerCase().contains("hugo")) {
 gp_id = "p" + id;
 if (!nw.nodes.containsKey(gp_id)) {
 i++;
 n2 = new Node(gp_id, "protein", "hugo", "10000008", n.ids);
 new_nodes.put(gp_id, n2);
 //nw.nodes.put(gp_id, n2);
 other_ids = new ArrayList();
 sourcedbentry = new ArrayList();
 new_int.put("man" + i, new Interaction("man" + i, other_ids, n, n2, "expression", "man", sourcedbentry, "manexpr", "-"));
 //nw.interactions.put("man" + i, new Interaction("man" + i, other_ids, n, n2, "expression", "man", sourcedbentry, "manexpr", "-"));
 } else {
 if (!nw.nodes.get(id).upnbrs.containsKey(nw.nodes.get(gp_id))) {
 i++;
 n2 = nw.nodes.get(gp_id);
 other_ids = new ArrayList();
 sourcedbentry = new ArrayList();
 // new_int.put("man" + i, new Interaction("man" + i, other_ids, n, n2, "expression", "man", sourcedbentry, "manexpr", "-"));
 }
 }
 }
 } 
 }


 public void add_missing_genes_and_products(Network nw) {
 System.out.println("+++++++++++add_missing_genes_and_products++++++++++");
 System.out.println("Before statistics: nodes " + nw.nodes.size() + " interactions " + nw.interactions.size());
 String id, type, id_type, tmp_id;
 Node n2;
 Interaction int1;
 int i = 0;
 List<String> other_ids;
 List<String[]> sourcedbentry;
 Map<String, Node> new_nodes = new HashMap();
 Map<String, Interaction> tmp_map;
 Map<String, Interaction> new_int = new HashMap();
 boolean found;
 Node n1;
 Node[] nn = new Node[2];
 List<String> lookedat = new ArrayList();
 for (Interaction tr : nw.interactions.values()) {
 nn[0] = tr.int1;
 nn[1] = tr.int2;
 for (Node n : nn) {
 found = false;
 id = n.id;
 type = n.type;
 id_type = n.id_type;
 if (type.toLowerCase().equals("protein") && id_type.toLowerCase().contains("hugo") && !lookedat.contains(id)) {
 tmp_id = id.substring(1);
 if (!nw.nodes.containsKey(tmp_id)) {
 i++;
 n2 = new Node(tmp_id, "gene", "hugo", "10000008", n.ids);
 other_ids = new ArrayList();
 sourcedbentry = new ArrayList();
 int1 = new Interaction("man" + i, other_ids, n2, n, "expression", "man", sourcedbentry, "manexpr", "-");
 n2.downnbrs.put(n.id, int1);
 nw.nodes.get(n.id).upnbrs.put(n2.id, int1);
 new_int.put("man" + i, int1);
 new_nodes.put(tmp_id, n2);
 //new_nodes.put(id)
 lookedat.add(id);
 lookedat.add(tmp_id);
 //new_int.put("man" + i, new Interaction("man" + i, other_ids, n2, n, "expression", "man", sourcedbentry, "manexpr", "-"));
 }
 if (nw.nodes.containsKey(tmp_id)) {
 // tmp_map = new HashMap();
 for (String ss : n.upnbrs.keySet()) {
 if (ss.equals(tmp_id)) {
 found = true;
 lookedat.add(id);
 lookedat.add(tmp_id);
 break;
 }
 }
 if (!found) {
 i++;
 other_ids = new ArrayList();
 sourcedbentry = new ArrayList();
 int1 = new Interaction("man" + i, other_ids, nw.nodes.get(tmp_id), n, "expression", "man", sourcedbentry, "manexpr", "-");
 // n.upnbrs.put(tmp_id, int1);
 // n2 = nw.nodes.get(tmp_id);
 // n2.downnbrs.put(n.id, int1);
 nw.nodes.get(tmp_id).downnbrs.put(n.id, int1);
 nw.nodes.get(n.id).upnbrs.put(tmp_id, int1);
 // tmp_map.put(tmp_id, int1);
 new_int.put("man" + i, int1);
 lookedat.add(id);
 lookedat.add(tmp_id);
 }
 //nw.nodes.get(n.id).upnbrs.putAll(tmp_map);
 //tmp_map = null;
 }
 }
 if (type.toLowerCase().equals("gene") && id_type.toLowerCase().contains("hugo") && !lookedat.contains(id)) {
 tmp_id = "p" + id;
 if (!nw.nodes.containsKey(tmp_id)) {
 i++;
 n2 = new Node(tmp_id, "protein", "hugo", "10000008", n.ids);
 other_ids = new ArrayList();
 sourcedbentry = new ArrayList();
 int1 = new Interaction("man" + i, other_ids, n, n2, "expression", "man", sourcedbentry, "manexpr", "-");
 n2.upnbrs.put(n.id, int1);
 nw.nodes.get(n.id).downnbrs.put(n2.id, int1);
 new_nodes.put(tmp_id, n2);
 new_int.put("man" + i, int1);
 lookedat.add(id);
 lookedat.add(tmp_id);
 }
 if (nw.nodes.containsKey(tmp_id)) {
 // tmp_map = new HashMap();
 for (String ss : n.downnbrs.keySet()) {
 if (ss.equals(tmp_id)) {
 found = true;
 lookedat.add(id);
 lookedat.add(tmp_id);
 break;
 }
 }
 if (!found) {
 i++;
 other_ids = new ArrayList();
 sourcedbentry = new ArrayList();
 int1 = new Interaction("man" + i, other_ids, n, nw.nodes.get(tmp_id), "expression", "man", sourcedbentry, "manexpr", "-");
 nw.nodes.get(tmp_id).upnbrs.put(n.id, int1);
 nw.nodes.get(n.id).downnbrs.put(tmp_id, int1);
 // tmp_map.put(tmp_id, int1);
 new_int.put("man" + i, int1);
 lookedat.add(id);
 lookedat.add(tmp_id);
 }
 // nw.nodes.get(n.id).downnbrs.putAll(tmp_map);
 // tmp_map = null;
 }
 }
 }
 }

 nw.nodes.putAll(new_nodes);
 nw.interactions.putAll(new_int);
 System.out.println("After statistics: nodes " + nw.nodes.size() + " interactions " + nw.interactions.size());
 System.out.println("Added nodes " + new_nodes.size() + " int  " + new_int.size());
 System.out.println("i " + i);
 }



 //      if (!visited_nodes.containsKey(n)){ //&& !n2.contains(n)) {
 path_type = return_node_type(network.nodes.get(n).type);
 //          try {
 if (visited_nodes.containsKey(n)) {
 tmp_path_type = visited_nodes.get(curn.current_seed).clone();
 } else {
 tmp_path_type = new int[]{0, 0, 0, 0};
 }
 //           } catch (NullPointerException e) {
 //               tmp_path_type = only_ppi.clone();
 //         }
 tmp_path_type[path_type] = 1;
 //|| !visited_nodes.containsKey(n)
 if (!Arrays.equals(only_ppi, tmp_path_type) || !tmp_visited_nodes.containsKey(n)) {
 tmp_visited_nodes.put(n, tmp_path_type);
 }


 */
/*

 public List<PathUnit> getLongPathDF(Network nw, String n1, List<String> n2, int max) {
 // System.out.println("++++++++++getLongPathDF+++++++++");
 List<String> path;
 List<String> seeds = new ArrayList();
 Stack<PathUnit> stack = new Stack();
 Stack<PathUnit> found = new Stack();
 PathUnit curn;
 if (nw.nodes.containsKey(n1)) {
 for (String n : nw.nodes.get(n1).downnbrs.keySet()) {
 path = new ArrayList();
 path.add(nw.nodes.get(n1).downnbrs.get(n).id);
 stack.push(new PathUnit(path, seeds, n, ""));
 }
 for (String n : nw.nodes.get(n1).revnbrs.keySet()) {
 path = new ArrayList();
 path.add(nw.nodes.get(n1).revnbrs.get(n).id);
 stack.push(new PathUnit(path, seeds, n, ""));
 }
 }

 while (!stack.empty()) {
 curn = stack.pop();
 if (curn.path.size() > max) {
 continue;
 }
 for (String end : n2) {
 if (curn.current_seed.equals(end)) {
 found.push(curn);
 //continue;
 }
 }
 for (String n : nw.nodes.get(curn.current_seed).revnbrs.keySet()) {
 if (!nw.nodes.get(curn.current_seed).revnbrs.get(n).id.equals(curn.path.get(curn.path.size() - 1))) {
 path = new ArrayList();
 path.addAll(curn.path);
 path.add(nw.nodes.get(curn.current_seed).revnbrs.get(n).id);
 stack.push(new PathUnit(path, seeds, n, ""));
 }
 }
 for (String n : nw.nodes.get(curn.current_seed).downnbrs.keySet()) {
 path = new ArrayList();
 path.addAll(curn.path);
 path.add(nw.nodes.get(curn.current_seed).downnbrs.get(n).id);
 stack.push(new PathUnit(path, seeds, n, ""));
 }

 }
 stack = null;
 path = null;
 System.gc();
 return found;
 }

 public List<PathUnit> getLongPathDF_linkedlist(Network network, List<String> n1, List<String> n2, int max, String prefix, BufferedWriter out) throws IOException {
 // System.out.println("++++++++++getLongPathDF+++++++++");
 PathUnit curn, tmp_curn;
 List<String> path, seeds, n2_wo_identical;
 List<PathUnit> tmp_punitlist;
 LinkedList<PathUnit> stack = new LinkedList(), found = new LinkedList();
 Map<String, Depth_List_of_PathUnit> tmp_visited_nodes, visited_nodes = new HashMap();
 Depth_List_of_PathUnit dpunit;
 int depth, current_hg = 0, id = 0;

 for (String beg : n1) {
 current_hg++;
 System.out.println("Current hit gene " + beg + "(" + current_hg + "/" + n1.size() + ")");
 n2_wo_identical = new ArrayList();
 for (String end : n2) {
 if (!beg.equals(end)) {
 n2_wo_identical.add(end);
 }
 }
 if (!visited_nodes.containsKey(beg)) {
 visited_nodes.put(beg, new Depth_List_of_PathUnit());
 }
 if (network.nodes.containsKey(beg)) {
 for (String n : network.nodes.get(beg).downnbrs.keySet()) {
 if (!visited_nodes.containsKey(n)) {
 visited_nodes.put(n, new Depth_List_of_PathUnit());
 }
 path = new ArrayList();
 path.add(network.nodes.get(beg).downnbrs.get(n).id);
 seeds = new ArrayList();
 seeds.add(n);
 stack.push(new PathUnit(path, seeds, n, ""));
 }
 for (String n : network.nodes.get(beg).revnbrs.keySet()) {
 if (!visited_nodes.containsKey(n)) {
 visited_nodes.put(n, new Depth_List_of_PathUnit());
 }
 path = new ArrayList();
 path.add(network.nodes.get(beg).revnbrs.get(n).id);
 seeds = new ArrayList();
 seeds.add(n);
 stack.push(new PathUnit(path, seeds, n, ""));
 }
 }

 while (!stack.isEmpty()) {
 curn = stack.poll();
 for (String end : n2_wo_identical) {
 if (curn.current_seed.equals(end)) {
 found.push(curn);
 for (int i = 0; i <= curn.seeds.size() - 1; i++) {
 tmp_punitlist = visited_nodes.get(curn.seeds.get(i)).punits;
 tmp_punitlist.add(curn);
 //tmp_punitlist.add(new PathUnit(curn.path.subList(i, curn.path.size() - 1), curn.seeds.subList(i, curn.seeds.size() - 1), curn.seeds.get(i)));
 depth = curn.path.size() - i - 1;
 dpunit = new Depth_List_of_PathUnit(depth, tmp_punitlist);
 visited_nodes.put(curn.seeds.get(i), dpunit); //continue;
 //visited_nodes.get(curn.seeds.get(i)).max_depth = curn.path.size() - i - 1;
 }
 }
 }

 if (curn.path.size() >= max) {
 for (int i = 0; i <= curn.seeds.size() - 1; i++) {
 visited_nodes.get(curn.seeds.get(i)).max_depth = curn.path.size() - i - 1;
 }
 continue;
 }

 for (String n : network.nodes.get(curn.current_seed).revnbrs.keySet()) {
 if (!network.nodes.get(curn.current_seed).revnbrs.get(n).id.equals(curn.path.get(curn.path.size() - 1))) {
 if (visited_nodes.containsKey(n) && (visited_nodes.get(n).max_depth + curn.path.size() >= max) && visited_nodes.get(n).punits.isEmpty()) {
 continue;
 }

 if (visited_nodes.containsKey(n) && !visited_nodes.get(n).punits.isEmpty() && (visited_nodes.get(n).max_depth + curn.path.size() >= max)) {
 /*     tmp_visited_nodes = new HashMap();
 for (PathUnit pu : visited_nodes.get(n).punits) {
 if (pu.path.size() + curn.path.size() <= max) {
 tmp_curn = new PathUnit();
 tmp_curn.path.addAll(curn.path);
 tmp_curn.path.addAll(pu.path);
 tmp_curn.seeds.addAll(curn.seeds);
 tmp_curn.seeds.addAll(pu.seeds);
 tmp_curn.current_seed = pu.current_seed;
 found.push(tmp_curn);
 for (int i = 0; i <= curn.seeds.size()-1; i++) {
 tmp_punitlist =  new ArrayList(); //
 tmp_punitlist.addAll(visited_nodes.get(tmp_curn.seeds.get(i)).punits);
 tmp_punitlist.add(new PathUnit(tmp_curn.path.subList(i , tmp_curn.path.size()-1), tmp_curn.seeds.subList(i , tmp_curn.seeds.size()-1), tmp_curn.seeds.get(i)));
 depth = visited_nodes.get(n).max_depth + curn.seeds.size() - i-1;
 dpunit = new Depth_List_of_PathUnit(depth, tmp_punitlist);
 tmp_visited_nodes.put(curn.seeds.get(i), dpunit);
 }
 }
 }
 visited_nodes.putAll(tmp_visited_nodes);
 continue; */
/*                       }

 path = new ArrayList();
 path.addAll(curn.path);
 path.add(network.nodes.get(curn.current_seed).revnbrs.get(n).id);
 seeds = new ArrayList();
 seeds.addAll(curn.seeds);
 seeds.add(n);
 stack.push(new PathUnit(path, seeds, n, ""));
 if (!visited_nodes.containsKey(n)) {
 visited_nodes.put(n, new Depth_List_of_PathUnit());
 }
 }
 }
 for (String n : network.nodes.get(curn.current_seed).downnbrs.keySet()) {
 if (visited_nodes.containsKey(n) && (visited_nodes.get(n).max_depth + curn.path.size() >= max) && visited_nodes.get(n).punits.isEmpty()) {
 continue;
 }
 if (visited_nodes.containsKey(n) && !visited_nodes.get(n).punits.isEmpty() && (visited_nodes.get(n).max_depth + curn.path.size() >= max)) {
 /*
 tmp_visited_nodes = new HashMap();
 for (PathUnit pu : visited_nodes.get(n).punits) {
 if (pu.path.size() + curn.path.size() <= max) {
 tmp_curn = new PathUnit();
 tmp_curn.path.addAll(curn.path);
 tmp_curn.path.addAll(pu.path);
 tmp_curn.seeds.addAll(curn.seeds);
 tmp_curn.seeds.addAll(pu.seeds);
 tmp_curn.current_seed = pu.current_seed;
 found.push(tmp_curn);
 for (int i = 0; i <= curn.seeds.size()-1; i++) {
 tmp_punitlist = new ArrayList();
 tmp_punitlist.addAll(visited_nodes.get(tmp_curn.seeds.get(i)).punits);
 tmp_punitlist.add(new PathUnit(tmp_curn.path.subList(i, tmp_curn.path.size()-1), tmp_curn.seeds.subList(i , tmp_curn.seeds.size()-1), tmp_curn.seeds.get(i)));
 depth = visited_nodes.get(n).max_depth + curn.seeds.size() - i-1;
 dpunit = new Depth_List_of_PathUnit(depth, tmp_punitlist);
 tmp_visited_nodes.put(curn.seeds.get(i), dpunit);
 }
 }
 }
 visited_nodes.putAll(tmp_visited_nodes);
 continue; */
/*                 }
 path = new ArrayList();
 path.addAll(curn.path);
 path.add(network.nodes.get(curn.current_seed).downnbrs.get(n).id);
 seeds = new ArrayList();
 seeds.addAll(curn.seeds);
 seeds.add(n);
 stack.push(new PathUnit(path, seeds, n, ""));
 if (!visited_nodes.containsKey(n)) {
 visited_nodes.put(n, new Depth_List_of_PathUnit());
 }
 }
 }

 for (int i = 0; i < found.size(); i++) {
 id++;
 out.write(prefix + id + "\t");
 out.write(beg);
 for (int j = 0; j < found.get(i).path.size(); j++) {
 out.write("\t" + found.get(i).path.get(j));
 }
 out.write("\t" + found.get(i).current_seed);
 out.newLine();
 }
 //out.flush();
 found = new LinkedList();

 }
 stack = null;
 path = null;

 System.gc();
 return found;
 }


 /*
 if (visited_nodes.containsKey(curn.current_seed)) {
 if (!visited_nodes.get(curn.current_seed).punits.isEmpty()) {
 for (PathUnit punit : visited_nodes.get(curn.current_seed).punits) {
 if ((curn.path.size() + punit.path.size()) <= max) {
 curn.path.addAll(punit.path);
 curn.current_seed = punit.current_seed;
 curn.seeds.addAll(punit.seeds);
 found.push(curn);
 }
 }
 } else {
 continue;
 }
 }
 */
/*   public void getLongPathsforListDF(Network nw, Map<String, String[]> hugo, int max, String outf, String prefix) throws IOException {
 System.out.println("++++++++++getPathforListDF+++++++++");
 BufferedWriter out = new BufferedWriter(new FileWriter(outf));
 List<NetworkUtils.PathUnit> paths;
 int id = 0;
 int current_hg = 0;
 // for (String s : hg.keySet()) {
 //current_hg++;
 //System.out.println("Current hit gene " + s + "(" + current_hg + "/" + hg.keySet().size() + ")");
 paths = getLongPathDF_linkedlist(nw, new ArrayList(hg.keySet()), new ArrayList(fpl.keySet()), max, prefix, out);
 /*for (int i = 0; i < paths.size(); i++) {
 id++;
 out.write(prefix + id + "\t");
 out.write("s");
 for (int j = 0; j < paths.get(i).path.size(); j++) {
 out.write("\t" + paths.get(i).path.get(j));
 }
 out.write("\t" + paths.get(i).current_seed);
 out.newLine();
 }
 out.flush();
 // } */
/*        out.close();
 }




 public Network mergeLight(Network hprd, Network tp, Network kegg, Network transmir, Network mirtarbase, Network tfacts) {
 System.out.println("\n +++++++++++++mergeLight++++++++++++ \n");
 Network merged = new Network();
 merged.nodes.putAll(hprd.nodes);
 merged.interactions.putAll(hprd.interactions);
 List<Map<String, Node>> nd_maps = new ArrayList();
 nd_maps.add(null);
 nd_maps.add(null);
 nd_maps.add(null);
 nd_maps.add(null);
 nd_maps.add(null);
 nd_maps.add(null);
 nd_maps.add(null);
 nd_maps.add(3, tp.nodes);
 nd_maps.add(4, kegg.nodes);
 nd_maps.add(5, transmir.nodes);
 nd_maps.add(6, mirtarbase.nodes);
 nd_maps.add(7, tfacts.nodes);

 List<Map<String, Interaction>> int_maps = new ArrayList();
 int_maps.add(null);
 int_maps.add(null);
 int_maps.add(null);
 int_maps.add(null);
 int_maps.add(null);
 int_maps.add(null);
 int_maps.add(null);
 int_maps.add(3, tp.interactions);
 int_maps.add(4, kegg.interactions);
 int_maps.add(5, transmir.interactions);
 int_maps.add(6, mirtarbase.interactions);
 int_maps.add(7, tfacts.interactions);

 int[] db = {3, 4, 5, 6, 7};
 Node n;
 //nodes
 for (int i : db) {
 for (String id : nd_maps.get(i).keySet()) {
 if (merged.nodes.containsKey(id)) {
 merged.nodes.get(id).db_flag = merged.nodes.get(id).db_flag + i;
 merged.nodes.get(id).id_type = merged.nodes.get(id).id_type + nd_maps.get(i).get(id).id_type;
 merged.nodes.get(id).ids.addAll(nd_maps.get(i).get(id).ids);
 for (String s : nd_maps.get(i).get(id).downnbrs.keySet()) {
 if (!merged.nodes.get(id).downnbrs.containsKey(s)) {
 merged.nodes.get(id).downnbrs.put(s, nd_maps.get(i).get(id).downnbrs.get(s));
 }
 }
 for (String s : nd_maps.get(i).get(id).upnbrs.keySet()) {
 if (!merged.nodes.get(id).upnbrs.containsKey(s)) {
 merged.nodes.get(id).upnbrs.put(s, nd_maps.get(i).get(id).upnbrs.get(s));
 }
 }
 for (String s : nd_maps.get(i).get(id).revnbrs.keySet()) {
 if (!merged.nodes.get(id).revnbrs.containsKey(s)) {
 merged.nodes.get(id).revnbrs.put(s, nd_maps.get(i).get(id).revnbrs.get(s));
 }
 }

 } else {
 merged.nodes.put(id, nd_maps.get(i).get(id));
 }
 }
 }
 //interactions
 boolean flag;
 Node db1;
 Node db2;
 Node m1;
 Node m2;
 String s1, s2, s3, s4;

 db = new int[]{3, 7};
 //db = new int[]{};  // чтобы быстрее считался при тестировании
 merged.interactions.putAll(kegg.interactions);
 merged.interactions.putAll(transmir.interactions);
 merged.interactions.putAll(mirtarbase.interactions);

 for (int i : db) {
 for (String id : int_maps.get(i).keySet()) {
 flag = false;
 db1 = int_maps.get(i).get(id).int1;
 db2 = int_maps.get(i).get(id).int2;
 for (String m_id : merged.interactions.keySet()) {
 m1 = merged.interactions.get(m_id).int1;
 m2 = merged.interactions.get(m_id).int2;
 if (db1.id.equals(m1.id) && db2.id.equals(m2.id)) {
 merged.interactions.get(m_id).other_ids.add(id);
 merged.interactions.get(m_id).sourcedb = merged.interactions.get(m_id).sourcedb + db_name_by_number(i);
 merged.interactions.get(m_id).sourcedbentry.add(int_maps.get(i).get(id).sourcedbentry.get(0));
 flag = true;
 break;
 //копия в любом случае
 }
 /* if (db1.id.equals(m2.id) && db2.id.equals(m1.id) && db1.type.toLowerCase().equals("protein") && db2.type.toLowerCase().equals("protein") && m1.type.toLowerCase().equals("protein") && m2.type.toLowerCase().equals("protein")) {
 merged.interactions.get(m_id).other_ids.add(id);
 merged.interactions.get(m_id).sourcedb = merged.interactions.get(m_id).sourcedb + db_name_by_number(i);
 merged.interactions.get(m_id).sourcedbentry.add(int_maps.get(i).get(id).sourcedbentry.get(0));
 flag = true;
 break;
 //перевернутый PPI   - копия
 } */
/*                 if (db1.type.equals("protein") && db2.type.equals("protein") && m1.type.equals("protein") && m2.type.equals("protein")) {
 }
 }
 if (!flag) {
 merged.interactions.put(id, int_maps.get(i).get(id));
 }
 }
 }
 System.out.println("Statistics:\n\n" + "Merged entities= " + merged.nodes.size() + " Merged interactions= " + merged.interactions.size());
 System.out.println("\nHPRD entities= " + hprd.nodes.size() + " HPRD interactions= " + hprd.interactions.size());
 System.out.println("TP entities= " + tp.nodes.size() + " TP interactions= " + tp.interactions.size());
 System.out.println("KEGG entities= " + kegg.nodes.size() + " KEGG interactions= " + kegg.interactions.size());
 System.out.println("transmir entities= " + transmir.nodes.size() + " transmir interactions= " + transmir.interactions.size());
 System.out.println("mirtarbase entities= " + mirtarbase.nodes.size() + " mirtarbase interactions= " + mirtarbase.interactions.size());
 System.out.println("tfacts entities= " + tfacts.nodes.size() + " tfacts interactions= " + tfacts.interactions.size());

 int nds = hprd.nodes.size() + tp.nodes.size() + kegg.nodes.size() + transmir.nodes.size() + mirtarbase.nodes.size() + tfacts.nodes.size();
 int ints = hprd.interactions.size() + tp.interactions.size() + kegg.interactions.size() + transmir.interactions.size() + mirtarbase.interactions.size() + tfacts.interactions.size();
 System.out.println("\nsumm= " + nds + "  interactions= " + ints + "\n");

 return merged;
 //end method
 }


 public Network merge_direct_indirect(Network n1, Network n2) {
 System.out.println("\n +++++++++++++merge direct indirect++++++++++++ \n");
 Network mrg = new Network();
 boolean match = false;
 int i = 0, j = 0, k = 0, l = 0;
 mrg.nodes.putAll(n2.nodes);
 mrg.nodes.putAll(n1.nodes);

 mrg.interactions.putAll(n1.interactions);
 for (Interaction n2int : n2.interactions.values()) {
 for (Interaction n1int : n1.interactions.values()) {
 if ((n1int.int1.id.equals(n2int.int1.id) && n1int.int2.id.equals(n2int.int2.id)) || (n1int.int1.id.equals(n2int.int2.id) && n1int.int2.id.equals(n2int.int1.id))) {
 match = true;
 i++;
 break;
 }
 }
 if (!match) {
 mrg.interactions.put(n2int.id, n2int);
 }
 match = false;
 }
 System.out.println("First network nodes " + n1.nodes.size() + " interactions " + n1.interactions.size());
 System.out.println("Second network nodes " + n2.nodes.size() + " interactions " + n2.interactions.size());
 System.out.println("Overlap interactions" + i);
 return mrg;
 }
 public void checkHPRD_to_network_utilities(Network all) throws FileNotFoundException, IOException {
 BufferedReader rd = new BufferedReader(new FileReader("F:\\Dropbox\\_Работа\\Programms\\DATA\\Muscle Differentiation\\miRNA results\\JAG1_all_paths.txt"));
 BufferedWriter wr2 = new BufferedWriter(new FileWriter("F:\\Dropbox\\_Работа\\Programms\\DATA\\Muscle Differentiation\\miRNA results\\JAG1_all_paths_detailed.txt"));
 String s, tmp, end;
 String[] ss;
 int[] arr;
 int count;
 //Map<String, String[]> freq = new HashMap();
 Map<String, EdgeAttributes> occur = new HashMap();
 // pHGNC:7611 , pHGNC:4223 , pHGNC:7566  , pHGNC:5466

 String[] fp = new String[]{"pHGNC:7611", "pHGNC:4223", "pHGNC:7566", "pHGNC:5466"};
 while ((s = rd.readLine()) != null) {
 ss = s.split("\t");
 tmp = ss[0] + "\t";
 for (int i = 1; i < ss.length - 1; i++) {
 if (ss[i].contains("hprd:")) {
 tmp = tmp + "\t" + all.interactions.get(ss[i]).int1.id + "-" + all.interactions.get(ss[i]).int2.id;
 }
 }
 wr2.write(s + "\t" + tmp);
 wr2.newLine();

 }
 rd.close();
 // wr.close();
 wr2.close();
 }

 public void getConnectivity(Network nw, String listf, String outf) throws FileNotFoundException, IOException {
 BufferedReader rd = new BufferedReader(new FileReader(listf));
 BufferedWriter wr = new BufferedWriter(new FileWriter(outf));
 String[] ss = new String[3];
 String s;
 int deg = 0;
 while ((s = rd.readLine()) != null) {
 ss = s.split("\t");
 deg = nw.nodes.get(ss[0]).downnbrs.size() + nw.nodes.get(ss[1]).upnbrs.size() + nw.nodes.get(ss[0]).revnbrs.size();
 wr.write(s + "\t" + deg);
 wr.newLine();
 }
 rd.close();
 wr.close();
 }
 */
