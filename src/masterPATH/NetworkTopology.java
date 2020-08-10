package masterpath.masterpath;


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

/**
 * NetworkTopology class contains methods to calculate network's topological properties
 * 
 * @author Natalia Rubanova
 */
public class NetworkTopology {

    /**
     *
     * @param putils
     * @param nutils
     * @param dbutils
     * @param all
     * @param musclef_com
     * @param fp
     * @param foundDF
     * @param foundDF_mirnas
     * @throws IOException
     */
    private void get_and_calculate_for_random(PathwayManager putils, NetworkManager nutils, DBManager dbutils, Network all, String musclef_com, String fp, String foundDF, String foundDF_mirnas) throws IOException {
        Random randomGenerator = new Random();
        Map<String, String> r_hl = new HashMap();
        int count = 1000;
        int hg_number = 82;
        List<String> all_hugo_ids = new ArrayList();
        for (String s : dbutils.hugo.keySet()) {
            all_hugo_ids.add(s);
        }
        int randomInt;
        BufferedWriter r_hl_file;
        for (int i = 0; i < count; i++) {
            System.out.println("# " + i);
            r_hl_file = new BufferedWriter(new FileWriter(musclef_com + "HitGenes.txt" + "_random"));
            for (int j = 0; j < hg_number; j++) {
                randomInt = randomGenerator.nextInt(dbutils.hprd_to_hugo.size());
                r_hl.put(all_hugo_ids.get(randomInt), "");
                r_hl_file.write(all_hugo_ids.get(randomInt) + "\n");
                System.out.print(all_hugo_ids.get(randomInt) + "\t");
            }
            r_hl_file.close();
            //   nutils.getLongPathsforListDF(all,  dbutils.hugo, 6, foundDF + "_random", "ququ");
            putils.find_miRNAs_on_pathways(foundDF + "_random", foundDF_mirnas + "_random", dbutils.mirtarbase, r_hl, nutils.fpl, 6, 2, "");
        }
    }




    private void hubStatistics(Network all, String outfile, int len) throws IOException {
        BufferedWriter wr = new BufferedWriter(new FileWriter(outfile));
        ArrayList<String> nodes = new ArrayList();
        ArrayList<String> nbs = new ArrayList();
        Random rand = new Random();
        Random rand1 = new Random();
        int nd, j = 20000, con;
        String node, next, scon;
        for (String s : all.nodes.keySet()) {
            nodes.add(s);
        }
        while (j > 0) {
            nd = rand.nextInt(all.nodes.size());
            node = nodes.get(nd);
            con = (all.nodes.get(node).downnbrs.size() + all.nodes.get(node).revnbrs.size()) * (all.nodes.get(node).upnbrs.size() + all.nodes.get(node).revnbrs.size());
            scon = (all.nodes.get(node).downnbrs.size() + all.nodes.get(node).revnbrs.size()) + "*" + (all.nodes.get(node).upnbrs.size() + all.nodes.get(node).revnbrs.size());
            next = node;
            for (int k = 0; k < len; k++) {
                for (String s : all.nodes.get(next).downnbrs.keySet()) {
                    nbs.add(s);
                }
                for (String s : all.nodes.get(next).revnbrs.keySet()) {
                    nbs.add(s);
                }
                nd = rand.nextInt(nbs.size());
                next = nbs.get(nd);
                if (con < (all.nodes.get(next).downnbrs.size() + all.nodes.get(next).revnbrs.size()) * (all.nodes.get(next).upnbrs.size() + all.nodes.get(next).revnbrs.size())) {
                    con = (all.nodes.get(next).downnbrs.size() + all.nodes.get(next).revnbrs.size()) * (all.nodes.get(next).upnbrs.size() + all.nodes.get(next).revnbrs.size());
                    scon = (all.nodes.get(next).downnbrs.size() + all.nodes.get(next).revnbrs.size()) + "*" + (all.nodes.get(next).upnbrs.size() + all.nodes.get(next).revnbrs.size());
                }
            }
            j--;
            wr.write(con + " " + scon + "\n");
            // System.out.println(j+ " " + con + " " + scon);
        }

        wr.close();
    }

    private Network remove_genes(Network all) {
        Network all_new = new Network();
        return all_new;
    }

    /**
     * Calculate average clustering coefficient for PPI network
     * @param all Network
     */
    public void calculate_average_clustering_coefficient_ppi(Network all) {
        double ccoef = 0.0;
        String n2 = "";
        double degree;
        double edges;
        double number = 0.0;
        double number_of_nbrs = 0.0;
        double total_number_of_edges = 0.0;
        double density = 0.0;
        double tmp;
        List<String> nbrs;
        for (Node n : all.nodes.values()) {
            if (n.type.equals("protein")) {
                if (n.revnbrs.isEmpty()) {
                    continue;
                }
                number++;
                degree = 0.0;
                edges = 0.0;
                degree = n.revnbrs.size();
                //System.out.println(degree);
                nbrs = new ArrayList();
                for (Interaction edg : n.revnbrs.values()) {
                    if (!edg.int1.id.equals(n.id)) {
                        n2 = edg.int1.id;
                    }
                    if (!edg.int2.id.equals(n.id)) {
                        n2 = edg.int2.id;
                    }
                    nbrs.add(n2);
                }
                number_of_nbrs = number_of_nbrs + nbrs.size();
                total_number_of_edges = total_number_of_edges + nbrs.size();

                for (String n3 : nbrs) {
                    for (Interaction edg : all.nodes.get(n3).revnbrs.values()) {
                        if (edg.int1.id.equals(n3) && nbrs.contains(edg.int2.id)) {
                            edges++;
                        }
                        if (edg.int2.id.equals(n3) && nbrs.contains(edg.int1.id)) {
                            edges++;
                        }
                    }
                }

                tmp = (2 * edges * 1.0) / (degree * (degree - 1));
                if (tmp > 0.0) {
                    ccoef = ccoef + tmp;
                }
            }
        }
        ccoef = ccoef / number;
        number_of_nbrs = number_of_nbrs / number;
        total_number_of_edges = total_number_of_edges / 2;
        density = (2 * total_number_of_edges) / (number * (number - 1));
        System.out.println("average clustering coefficient = " + ccoef);
        System.out.println("average number of neighbours = " + number_of_nbrs);
        System.out.println("density = " + density);
    }

      /**
     * Calculate average clustering coefficient for direct network
     * @param all Network
     */
    public void calculate_average_clustering_coefficient_direct(Network all) {
        double ccoef = 0.0;
        List<String> nodes = new ArrayList();
        for (Node n : all.nodes.values()) {
            if ((!n.revnbrs.isEmpty() || !n.upnbrs.isEmpty() || !n.downnbrs.isEmpty()) && !n.type.equals("gene")) {
                nodes.add(n.id);
            }
        }
        byte[][] matrix = new byte[nodes.size()][nodes.size()];
        for (int i = 0; i < nodes.size(); i++) {
            for (int j = 0; j < nodes.size(); j++) {
                if (i == j) {
                    matrix[i][j] = 0;
                } else {
                    matrix[i][j] = 0;
                }
            }
        }
        int f, s;
        int test = 0;
        int number_of_edges = 0;
        int number_of_edges_direct = 0;
        int number_of_edges_indirect = 0;
        int number_of_genes = 0;
        List<String> visited_genes = new ArrayList();
        for (Interaction edg : all.interactions.values()) {
            if (!edg.int1.id.equals(edg.int2.id)) {
                if (edg.int1.type.toLowerCase().equals("gene")) {
                    f = nodes.indexOf("p" + edg.int1.id);
                    if (!visited_genes.contains(edg.int1.id)) {
                        number_of_genes++;
                        visited_genes.add(edg.int1.id);
                    }
                } else {
                    f = nodes.indexOf(edg.int1.id);
                }
                if (edg.int2.type.toLowerCase().equals("gene")) {
                    s = nodes.indexOf("p" + edg.int2.id);
                    if (!visited_genes.contains(edg.int2.id)) {
                        number_of_genes++;
                        visited_genes.add(edg.int2.id);
                    }
                } else {
                    s = nodes.indexOf(edg.int2.id);
                }

                if (f > 0 && s > 0 && f != s && matrix[f][s] == 0) {
                    if (edg.id.startsWith("man")) {
                        System.out.println("ACTUNG");
                        break;
                    }
                    if (edg.dir.equalsIgnoreCase("indirect")) {
                        matrix[f][s] = 1;
                        matrix[s][f] = 1;
                        number_of_edges = number_of_edges + 2;
                        number_of_edges_indirect++;
                    } else {
                        matrix[f][s] = 1;
                        number_of_edges++;
                        number_of_edges_direct++;
                    }
                }
            }
        }
        int number = 0;
        double tmp_ccoef;
        int d_rev;
        int d_total;
        int network_degree = 0;
        System.out.println(" number of genes " + number_of_genes + " total nodes " + nodes.size());
        System.out.println(" number of edges " + number_of_edges);
        System.out.println(" number of edges direct" + number_of_edges_direct);
        System.out.println(" number of edges indirect" + number_of_edges_indirect);
        List<Integer> j_list = new ArrayList();
        for (int i = 0; i < nodes.size(); i++) {
            //System.out.println(i);
            number++;
            tmp_ccoef = 0.0;
            d_total = 0;
            d_rev = 0;

            for (int j = 0; j < nodes.size(); j++) {
                if (i == j || (matrix[i][j] + matrix[j][i] == 0)) {
                    continue;
                }
                d_total = d_total + matrix[i][j] + matrix[j][i];
                network_degree = network_degree + matrix[i][j] + matrix[j][i];
                d_rev = d_rev + matrix[i][j] * matrix[j][i];
                for (int k = 0; k < nodes.size(); k++) {
                    //System.out.println(i + " " + j + " " + k);
                    if (k == j) {
                        continue;
                    }
                    tmp_ccoef = tmp_ccoef + (matrix[i][j] + matrix[j][i]) * (matrix[i][k] + matrix[k][i]) * (matrix[k][j] + matrix[j][k]);
                }
            }
            if (d_total * (d_total - 1) - 2 * d_rev != 0) {
                tmp_ccoef = tmp_ccoef / (2 * (d_total * (d_total - 1) - 2 * d_rev));
            } else {
                tmp_ccoef = 0;
            }
            ccoef = ccoef + tmp_ccoef;
            //System.out.println(ccoef);
        }
        ccoef = ccoef / number;
        double density = (1.0 * number_of_edges) / (nodes.size() * (nodes.size() - 1));
        double number_of_nbrs = (1.0 * network_degree) / (2 * nodes.size());
        System.out.println("average clustering coefficient = " + ccoef);
        System.out.println("average number of neighbours = " + number_of_nbrs);
        System.out.println("density = " + density);
    }
    
    
    /**
     * Calculate number of connected components PPI network 
     * @param all Network
     */
    public void calculate_number_of_connected_components_ppi(Network all) {
        List<List<String>> components = new ArrayList();
        int number_of_components = 0;
        List<String> nbrs;
        String n2 = "";
        int path = 0;
        boolean add = false;
        List<String> new_level;
        List<String> visited_nodes = new ArrayList();
        for (Node n : all.nodes.values()) {
            if (n.type.equals("protein") && !visited_nodes.contains(n.id)) {
                if (n.revnbrs.isEmpty()) {
                    continue;
                }
                path = 0;
                number_of_components++;
                nbrs = new ArrayList();
                visited_nodes.add(n.id);
                for (Interaction edg : n.revnbrs.values()) {
                    if (!edg.int1.id.equals(n.id) && !visited_nodes.contains(edg.int1.id)) {
                        n2 = edg.int1.id;
                        nbrs.add(n2);
                        visited_nodes.add(n2);
                    }
                    if (!edg.int2.id.equals(n.id) && !visited_nodes.contains(edg.int2.id)) {
                        n2 = edg.int2.id;
                        nbrs.add(n2);
                        visited_nodes.add(n2);
                    }
                }
                new_level = new ArrayList();
                new_level.addAll(nbrs);

                if (nbrs.isEmpty()) {
                    // System.out.println("A " + n.id + " " + n.revnbrs.size() );
                }
                while (!new_level.isEmpty()) {
                    path++;
                    for (String n3 : nbrs) {
                        new_level.remove(n3);
                        for (Interaction edg : all.nodes.get(n3).revnbrs.values()) {
                            if (!visited_nodes.contains(edg.int1.id)) {
                                new_level.add(edg.int1.id);
                                visited_nodes.add(edg.int1.id);
                            }
                            if (!visited_nodes.contains(edg.int2.id)) {
                                new_level.add(edg.int2.id);
                                visited_nodes.add(edg.int2.id);
                            }
                        }
                    }
                    nbrs.clear();
                    nbrs.addAll(new_level);
                }
                System.out.println("path = " + path);
            }
        }
        System.out.println("number of connected components = " + number_of_components);
    }

    /**
     * Calculate number of connected components for direct network
     * @param all Network
     */
    public void calculate_number_of_connected_components_direct(Network all) {
        int number_of_components = 0;
        List<String> nbrs;
        String n2;
        int path = 0;
        List<String> new_level;
        List<String> visited_nodes = new ArrayList();
        List<Interaction> all_nbrs = new ArrayList();
        for (Node n : all.nodes.values()) {
            if (!visited_nodes.contains(n.id)) {
                if (n.revnbrs.isEmpty() && n.upnbrs.isEmpty() && n.downnbrs.isEmpty()) {
                    continue;
                }
                path = 0;
                number_of_components++;
                nbrs = new ArrayList();
                visited_nodes.add(n.id);
                all_nbrs.clear();
                all_nbrs.addAll(n.revnbrs.values());
                all_nbrs.addAll(n.upnbrs.values());
                all_nbrs.addAll(n.downnbrs.values());
                for (Interaction edg : all_nbrs) {
                    if (!edg.int1.id.equals(n.id) && !visited_nodes.contains(edg.int1.id)) {
                        n2 = edg.int1.id;
                        nbrs.add(n2);
                        visited_nodes.add(n2);
                    }
                    if (!edg.int2.id.equals(n.id) && !visited_nodes.contains(edg.int2.id)) {
                        n2 = edg.int2.id;
                        nbrs.add(n2);
                        visited_nodes.add(n2);
                    }
                }
                new_level = new ArrayList();
                new_level.addAll(nbrs);

                if (nbrs.isEmpty()) {
                    // System.out.println("A " + n.id + " " + n.revnbrs.size() );
                }

                while (!new_level.isEmpty()) {
                    path++;
                    for (String n3 : nbrs) {
                        new_level.remove(n3);

                        all_nbrs.clear();
                        all_nbrs.addAll(all.nodes.get(n3).revnbrs.values());
                        all_nbrs.addAll(all.nodes.get(n3).upnbrs.values());
                        all_nbrs.addAll(all.nodes.get(n3).downnbrs.values());
                        for (Interaction edg : all_nbrs) {
                            if (!visited_nodes.contains(edg.int1.id)) {
                                new_level.add(edg.int1.id);
                                visited_nodes.add(edg.int1.id);
                            }
                            if (!visited_nodes.contains(edg.int2.id)) {
                                new_level.add(edg.int2.id);
                                visited_nodes.add(edg.int2.id);
                            }
                        }
                    }
                    nbrs.clear();
                    nbrs.addAll(new_level);
                }
                //System.out.println("path = " + path);
            }
        }
        System.out.println("number of connected components = " + number_of_components);
    }
    
   
    /**
     * Calculate diameter
     * @param all Network
     */
    public void calculate_diameter_ppi(Network all) {
        List<String> nodes = new ArrayList();
        for (Node n : all.nodes.values()) {
            if (n.type.equals("protein") && !n.revnbrs.isEmpty()) {
                nodes.add(n.id);
            }
        }
        int[][] matrix = new int[nodes.size()][nodes.size()];
        for (int i = 0; i < nodes.size(); i++) {
            for (int j = 0; j < nodes.size(); j++) {
                if (i == j) {
                    matrix[i][j] = 0;
                } else {
                    matrix[i][j] = 10000;
                    //System.out.println( matrix[i][j]);
                }
            }
        }
        int f, s;
        for (Interaction edg : all.interactions.values()) {
            if (edg.int1.type.equals("protein") && edg.int2.type.equals("protein") && !edg.int1.id.equals(edg.int2.id)) {
                f = nodes.indexOf(edg.int1.id);
                s = nodes.indexOf(edg.int2.id);
                matrix[f][s] = 1;
                matrix[s][f] = 1;
            }
        }
        for (int k = 0; k < nodes.size(); k++) {
            for (int i = 0; i < nodes.size(); i++) {
                for (int j = 0; j < nodes.size(); j++) {
                    if (matrix[i][j] > matrix[i][k] + matrix[k][j]) {
                        matrix[i][j] = matrix[i][k] + matrix[k][j];
                        //System.out.println(matrix[i][j]);
                    }
                }
            }
        }

        int diametr = 0;
        for (int i = 0; i < nodes.size(); i++) {
            for (int j = 0; j < nodes.size(); j++) {
                if (diametr < matrix[i][j] && matrix[i][j] < 10000) {
                    diametr = matrix[i][j];
                }
            }
        }
        int aspl_number = 0;
        double aspl = 0.0;
        for (int i = 0; i < nodes.size(); i++) {
            for (int j = 0; j < nodes.size(); j++) {
                if (i != j && matrix[i][j] < 1000) {
                    aspl = aspl + matrix[i][j];
                    aspl_number++;
                }
            }
        }
        aspl = aspl / aspl_number;
        System.out.println("diametr = " + diametr);
        System.out.println("aspl = " + aspl);
    }

    /**
     * Calculate diameter of a direct network 
     * @param all Network
     */
    public void calculate_diameter_direct(Network all) {
        List<String> nodes = new ArrayList();
        for (Node n : all.nodes.values()) {
            if ((!n.revnbrs.isEmpty() || !n.upnbrs.isEmpty() || !n.downnbrs.isEmpty()) && !n.type.equals("gene")) {
                nodes.add(n.id);
            }
        }
        System.out.println(" nodes size " + nodes.size());
        byte[][] matrix = new byte[nodes.size()][nodes.size()];
        for (int i = 0; i < nodes.size(); i++) {
            for (int j = 0; j < nodes.size(); j++) {
                if (i == j) {
                    matrix[i][j] = 0;
                } else {
                    matrix[i][j] = 60;
                }
            }
        }
        int f, s;
        for (Interaction edg : all.interactions.values()) {
            if (!edg.int1.id.equals(edg.int2.id)) {
                if (edg.int1.type.toLowerCase().equals("gene")) {
                    f = nodes.indexOf("p" + edg.int1.id);
                } else {
                    f = nodes.indexOf(edg.int1.id);
                }
                if (edg.int2.type.toLowerCase().equals("gene")) {
                    s = nodes.indexOf("p" + edg.int2.id);
                } else {
                    s = nodes.indexOf(edg.int2.id);
                }

                if (f > 0 && s > 0) {
                    if (edg.dir.equalsIgnoreCase("indirect")) {
                        matrix[f][s] = 1;
                        matrix[s][f] = 1;

                    } else {
                        matrix[f][s] = 1;

                    }
                }
            }
        }
        for (int k = 0; k < nodes.size(); k++) {
            System.out.println(k);
            for (int i = 0; i < nodes.size(); i++) {
                if (matrix[i][k] == 60) {
                    continue;
                }
                for (int j = 0; j < nodes.size(); j++) {
                    if (matrix[i][j] > matrix[i][k] + matrix[k][j]) {
                        matrix[i][j] = (byte) (matrix[i][k] + matrix[k][j]);
                    }
                }
            }
        }

        int diametr = 0;
        for (int i = 0; i < nodes.size(); i++) {
            for (int j = 0; j < nodes.size(); j++) {
                if (diametr < matrix[i][j] && matrix[i][j] < 60) {
                    diametr = matrix[i][j];
                }
            }
        }
        int aspl_number = 0;
        double aspl = 0.0;
        for (int i = 0; i < nodes.size(); i++) {
            for (int j = 0; j < nodes.size(); j++) {
                if (i != j && matrix[i][j] < 60) {
                    aspl = aspl + matrix[i][j];
                    aspl_number++;
                }
            }
        }
        aspl = aspl / aspl_number;
        System.out.println("diametr = " + diametr);
        System.out.println("aspl = " + aspl);
    }
    
    
    /**
     * Build degree distribution for PPI network
     * @param all Network
     * @param outfile Output file
     */
    public void get_degree_distribution_ppi(Network all, String outfile) throws IOException {
        BufferedWriter wr = new BufferedWriter(new FileWriter(outfile));
        int degree = 0;
        for (Node n : all.nodes.values()) {
            if (n.type.equals("protein")) {
                if (n.revnbrs.isEmpty()) {
                    continue;
                }
                System.out.println("" + n.revnbrs.size());
                wr.write("" + n.revnbrs.size());
                wr.newLine();
            }
        }
        wr.close();
    }

    /**
     * Build degree distribution for direct network
     * @param all Network
     * @param outfile Output file
     */
    public void get_degree_distribution_direct(Network all, String outfile) throws IOException {
        BufferedWriter wr_in = new BufferedWriter(new FileWriter(outfile + "_in"));
        BufferedWriter wr_out = new BufferedWriter(new FileWriter(outfile + "_out"));
        BufferedWriter wr_total = new BufferedWriter(new FileWriter(outfile + "_total"));
        int degree_in = 0, degree_out = 0, degree_total = 0;
        List<String> visited_nodes = new ArrayList();
        String gene_name, protein_name;
        for (Node n : all.nodes.values()) {
            if (!n.revnbrs.isEmpty() || !n.upnbrs.isEmpty() || !n.downnbrs.isEmpty()) {
                if (visited_nodes.contains(n.id)) {
                    continue;
                }
                if (n.type.equals("gene")) {
                    try {
                        gene_name = n.id;
                        protein_name = all.nodes.get("p" + n.id).id;
                    } catch (NullPointerException e) {
                        continue;
                    }
                    degree_in = n.upnbrs.size() + all.nodes.get("p" + n.id).upnbrs.size() + all.nodes.get("p" + n.id).revnbrs.size();
                    wr_in.write("" + degree_in + "\n");
                    degree_out = all.nodes.get("p" + n.id).downnbrs.size() + all.nodes.get("p" + n.id).revnbrs.size();
                    wr_out.write("" + degree_out + "\n");
                    degree_total = degree_in + degree_out;
                    wr_total.write("" + degree_total + "\n");
                    visited_nodes.add(n.id);
                    visited_nodes.add("p" + n.id);
                    continue;
                }
                if (n.type.equals("protein")) {
                    try {
                        gene_name = all.nodes.get(n.id.substring(1)).id;
                        protein_name = n.id;
                    } catch (NullPointerException e) {
                        continue;
                    }
                    degree_in = n.upnbrs.size() + n.revnbrs.size() + all.nodes.get(n.id.substring(1)).upnbrs.size();
                    wr_in.write("" + degree_in + "\n");
                    degree_out = n.downnbrs.size() + n.revnbrs.size();
                    wr_out.write("" + degree_out + "\n");
                    degree_total = degree_in + degree_out;
                    wr_total.write("" + degree_total + "\n");
                    visited_nodes.add(n.id);
                    visited_nodes.add(n.id.substring(1));
                    continue;
                }
                degree_in = n.revnbrs.size() + n.upnbrs.size();
                wr_in.write("" + degree_in + "\n");

                degree_out = n.revnbrs.size() + n.downnbrs.size();
                wr_out.write("" + degree_out + "\n");

                degree_total = 2 * n.revnbrs.size() + n.upnbrs.size() + n.downnbrs.size();
                wr_total.write("" + degree_total + "\n");

            }
        }
        wr_in.close();
        wr_out.close();
        wr_total.close();
    }
    
    
  
    

    
    

    /**
     * Calculate subnetwork properties
     * @param all Network 
     * @param fname File with pathways to create subnetwork
     */
    public void get_subnetwork_statistics(Network all, String fname) throws FileNotFoundException, IOException {
        BufferedReader rd = new BufferedReader(new FileReader(fname));
        String s;
        String[] ss;
        List<String> edges = new ArrayList();
        List<String> nodes = new ArrayList();
        List<String> genes = new ArrayList();
        int min_length = 10;
        int max_length = 0;
        while ((s = rd.readLine()) != null) {
            ss = s.split("\t");
            if (ss.length - 3 < min_length) {
                min_length = ss.length - 3;
            }
            if (ss.length - 3 > max_length) {
                max_length = ss.length - 3;
            }
            for (int i = 2; i < ss.length - 1; i++) {
                if (!edges.contains(ss[i])) {
                    edges.add(ss[i]);
                }
                if (!nodes.contains(all.interactions.get(ss[i]).int1.id)) {
                    nodes.add(all.interactions.get(ss[i]).int1.id);
                }
                if (!nodes.contains(all.interactions.get(ss[i]).int2.id)) {
                    nodes.add(all.interactions.get(ss[i]).int2.id);
                }
                if (all.interactions.get(ss[i]).int1.type.toLowerCase().equals("gene") && !genes.contains(all.interactions.get(ss[i]).int1.id)) {
                    genes.add(all.interactions.get(ss[i]).int1.id);
                }
                if (all.interactions.get(ss[i]).int2.type.toLowerCase().equals("gene") && !genes.contains(all.interactions.get(ss[i]).int2.id)) {
                    genes.add(all.interactions.get(ss[i]).int2.id);
                }
            }
        }
        System.out.println(" edges " + edges.size() + " nodes " + nodes.size() + " gene " + genes.size());
        System.out.println(" min length " + min_length + " max_length " + max_length);
    }
}
