package masterPATH;

import java.util.List;

/**
 * Interaction class describes an Interaction object
 * 
 * @author Natalia Rubanova
 */
public class Interaction {

    String id;
    List<String> other_ids;
    Node int1;
    Node int2;
    String type;
    String sourcedb;
    List<String> sourcedbentry;
    String quality;
    String dir;

    /**
     * Constructor 
     * @param id id of the interaction
     * @param other_ids ids from the source databases
     * @param int1 first interactor
     * @param int2 second interactor
     * @param type type of the interaction
     * @param sourcedb source database
     * @param sourcedbentry entry in the source database
     * @param quality confidence parameter of the interaction
     * @param dir direct or indirect interaction
     */
    public Interaction(String id, List<String> other_ids, Node int1, Node int2, String type, String sourcedb, List<String> sourcedbentry, String quality, String dir) {
        this.id = id;
        this.other_ids = other_ids;
        this.int1 = int1;
        this.int2 = int2;
        this.quality = quality;
        this.type = type;
        this.sourcedb = sourcedb;
        this.sourcedbentry = sourcedbentry;
        this.dir = dir;
    }

    @Override
    public String toString() {
        String s;
        s = id + "\t";
        for (String t : this.other_ids) {
            s = s + t + ";";
        }
        s = s + "\t" + int1.id + "\t" + int2.id + "\t" + type + "\t" + sourcedb + "\t";
        for (String t : this.sourcedbentry) {
            t = t.replaceAll("\t", "_");
            s = s + t + ";";
        }
        s = s + "\t" + quality + "\t" + dir;
        return s;
    }

    /**
     *Print all information about interaction in one line
     */
    public void print() {
        System.out.println(id + "\t" + int1.id + "\t" + int2.id + "\t" + type + "\t" + sourcedb + "\t" + quality + "\t" + dir);
    }
}
