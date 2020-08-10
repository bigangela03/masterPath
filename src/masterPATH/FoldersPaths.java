package masterpath.masterpath;

/**
 * FolderPaths class contains system paths to databases files, files with pathways, paths, nodes, etc
 * 
 * @author Natalia Rubanova
 */
public class FoldersPaths {

   
    static final String dbpath = "";
    static final String syspath = "";
    static final String sep = "";
     


    static final String hprd_folder = "" + sep;
    static final String transpath_folder = "" + sep;
    static final String kegg_folder = "" + sep;
    static final String transmir_folder = "" + sep;
    static final String mirtarbae_folder = "" + sep;
    static final String tfacts_folder = "" + sep;
    static final String hugo_folder = "" + sep;
    static String[] hugo_path = {dbpath + hugo_folder + "hugo_all_17Jun.txt"};

    static final String hp_folder = "" + sep;
    static final String sg_folder = "" + sep;
    static final String sl_folder = "" + sep;

    static String hugo_with_entrez = dbpath + hugo_folder + "HUGO_with_entrezid.txt";

    static final String hp_int = dbpath + hp_folder + "hippie_current.txt";
    static final String hp_int_high = dbpath + hp_folder + "hippie_high.txt";
    static final String hp_conv = dbpath + hp_folder + "hippie_hgnc_to_entrz.txt";

    static final String sg_int = dbpath + sg_folder + "all_data_edt.tsv";
    static final String sg_conv_uniprot_hgnc = dbpath + sg_folder + "signor_uniprot_hgnc.txt";

    static final String sl_int = dbpath + sl_folder + "06162016-signalink-VA5Eui_edt.csv";
    static final String sl_conv_uniprot_hgnc = dbpath + sl_folder + "signalink_uniprot_hgnc.txt";

    static String hprd_ent = dbpath + hprd_folder + "PROTEIN_NOMENCLATURE_EDT.txt";
    static String hprd_int = dbpath + hprd_folder + "BINARY_PROTEIN_PROTEIN_INTERACTIONS_EDT.txt";
    static String hprd_bioDB = dbpath + hprd_folder + "missed_hprd_to_hugo_bioDB.txt";
    static String hprd_conv = dbpath + hprd_folder + "hprd_id_conv.txt";


    static String kegg_ent_comp = dbpath + kegg_folder + "list_of_all_compounds.txt";
    static String kegg_ent_glyc = dbpath + kegg_folder + "list_of_all_glycans.txt";
    static String kegg_int = dbpath + kegg_folder + "mergedKEGGinteractions.txt";
    static String trsmir = dbpath + transmir_folder + "transmir_v1.2.txt";
    static String trsmir_crstbl = dbpath + transmir_folder + "miRNA_hugo.txt";
    static String mirtarb = dbpath + mirtarbae_folder + "hsa_MTI.txt";
    static String entrz_wo_hugo = dbpath + mirtarbae_folder + "entrez_wo_hugo.txt";
    static String tfactsf = dbpath + tfacts_folder + "tfacts.txt";
    static String tfacts_crosstable = dbpath + tfacts_folder + "tfacts_crosstable.txt";



    

    static final String stat = "Statistics" + sep;
    static final String inputstatf = "input" + sep;
    static final String outputstatf = "output" + sep;
    static final String outputstat = syspath + stat + outputstatf;

}


   