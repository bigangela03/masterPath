package masterPATH;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 *
 * @author Natalia Rubanova
 */
public class Node {

    String id;  //hugo id
    String type;
    String id_type;
    String db_flag;
    Map<String, Interaction> upnbrs;
    Map<String, Interaction> downnbrs;
    Map<String, Interaction> revnbrs;
    List<String[]> ids;

    /**
     *
     * @param id id of the Node
     * @param type type of the Node
     * @param id_type type of the id(nomenclature)
     * @param db_flag flag of the source database
     * @param ids ids from other nomenclatures
     */
    public Node(String id, String type, String id_type, String db_flag, List<String[]> ids) {
        this.id = id;
        this.upnbrs = new HashMap();
        this.downnbrs = new HashMap();
        this.revnbrs = new HashMap();
        this.ids = ids;
        this.type = type;
        this.db_flag = db_flag;
        this.id_type = id_type;
    }    
    
    @Override
    public String toString() {
        String s;
        s = id + "\t" +type+ "\t" +id_type+ "\t" +db_flag+"\t";
        for (String[] t : this.ids) {            
            for (String tt:t){
                s=s+"_"+tt;
            }
            s = s + ";";
        }
        return s;
    }
    
    /**
     * Print all node information in one line
     */
    public void print() {
        System.out.println(id + "\t" + type + "\t" + id_type + "\t" + db_flag + "\t" + ids.get(0));
    }
}
