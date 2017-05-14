package winbasecluster;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import edu.uci.ics.jung.graph.DirectedGraph;
import edu.uci.ics.jung.graph.DirectedSparseGraph;
import java.util.Map;
import java.util.Scanner;

public class Winbasecluster {
    
    static final String filepath = "clusterNumber10_step15\\cluster";
    static final String filepathold = "windowBasedClusterin\\cluster";
   static final int genenum = 182;
    static final int startwindow = 49;
    static final int endwindow = 97;
    static final String matrix = "matrix.txt";
    static final String group = "group.txt";
    static final String table = "table.txt";
    static final String connect = "connect.txt";
    static final String connectnew = "connectnew.txt";
    static final String xy = "xy.txt";
    static final String node = "node.txt";
    static HashMap<Integer, gene_group> genegroupID = new HashMap<>();
    static DirectedGraph<Integer, String> genegroup_trees = new DirectedSparseGraph<>();
    static HashMap<String, HashMap<Integer, List<Integer>>> genemap = new HashMap<>();
    static HashMap<Integer, List<Integer>> children = new HashMap<>();
    static HashMap<Integer, List<Integer>> parent = new HashMap<>();
    
    static BufferedReader readFile(String filename) throws FileNotFoundException {
        FileReader fr = new FileReader(filename);
        BufferedReader br = new BufferedReader(fr);
        return br;
    }
    
    public static void main(String[] args) throws FileNotFoundException, IOException {

        //readfile to int[][] gene_cluster
        FileOutputStream fosm = new FileOutputStream(new File(matrix));
        FileOutputStream fost = new FileOutputStream(new File(table));
        FileOutputStream fosc = new FileOutputStream(new File(connect));
        FileOutputStream fosg = new FileOutputStream(new File(group));
         FileOutputStream foscn = new FileOutputStream(new File(connectnew));
         FileOutputStream fosxy = new FileOutputStream(new File(xy));
        FileOutputStream fosnode = new FileOutputStream(new File(node));
        
        HashMap<Integer, gene_group> genegroupIDold = new HashMap<>();
              
        /*Scanner sc = new Scanner(System.in);
        System.out.println("filename: ");
        final String filepath = sc.nextLine();
        System.out.println("genenum: ");
        final int genenum = Integer.valueOf(sc.nextLine());
        System.out.println("startwindow: ");
        final int startwindow = Integer.valueOf(sc.nextLine());
        System.out.println("endwindow: ");
        final int endwindow = Integer.valueOf(sc.nextLine());*/
        int[][] gene_cluster = new int[genenum][endwindow - startwindow + 1];

        for (int i = 0; i <= endwindow - startwindow; i++) {
            int j = i + startwindow;
            /*if (j == 55 || j == 91 || j == 92 || j == 94) {
                for (int x = 0; x < genenum; x++) {
                    gene_cluster[x][i] = gene_cluster[x][i - 1];
                }
                continue;
            }*/
            String filepath0 = filepath + j + ".txt";
            BufferedReader file = readFile(filepath0);
            String line;
            
            int k = 0;
            while ((line = file.readLine()) != null) {              
                if (!line.startsWith("AT")) {
                    continue;
                }
                String[] lines = line.split("\t");
                gene_cluster[k++][i] = Integer.parseInt(lines[1].trim());
            }
        }
        
        //write file gene_cluster to matrix
        for (int i = 0; i < genenum; i++) {
            fosm.write((i + ":\t").getBytes());
            for (int j = 0; j < endwindow - startwindow; j++) {
                fosm.write((gene_cluster[i][j] + " ").getBytes());
            }
            fosm.write(("\n").getBytes());
        }
        System.out.println("read and write file succeed");

        

        //group and combination
        for(int sw = 0; sw < endwindow - startwindow; sw++){  
            HashMap<Integer, List<Integer>>genegrouplist = new HashMap<>();                        
            double used[] = new double[genenum];
            for(int i = 0; i < genenum; i++)used[i] = 0;         
            for(int genei = 0; genei < genenum; genei++){
                if(used[genei] == 1)continue;
                List<Integer>genelist = new ArrayList<>();
                
                for(int genej = genei+1; genej < genenum; genej++){
                    
                    if(gene_cluster[genei][sw] == gene_cluster[genej][sw] && gene_cluster[genei][sw+1] == gene_cluster[genej][sw+1]){
                       genelist.add(genej);
                       used[genej] = 1;
                    }     
                }
                if(genelist.isEmpty())continue;
                genegrouplist.put(genei, genelist);
            }           
            if(genegrouplist.isEmpty())continue;
            genemap.put(sw + "_" + (sw+1), genegrouplist);            
            
            for(int ew = sw+2; ew <= endwindow- startwindow; ew++){
                HashMap<Integer, List<Integer>>genegrouplist0 = genemap.get(sw + "_" + (ew-1));
                HashMap<Integer, List<Integer>>genegrouplist1 = new HashMap<>();
                
                for(int genei = 0; genei < genenum; genei++){
                    List<Integer> genelist0 = genegrouplist0.get(genei);
                    List<Integer> genelist1 = new ArrayList<>();
                    List<Integer> genelist2 = new ArrayList<>();
                    if(genelist0 == null || genelist0.isEmpty())continue;
                    for(Iterator<Integer> it0 = genelist0.iterator(); it0.hasNext();){
                        int gene0 = it0.next();
                        if(gene_cluster[genei][ew] == gene_cluster[gene0][ew]){
                            genelist1.add(gene0);
                            genelist2.add(gene0);
                        }
                        
                    }    
                    if(!genelist1.isEmpty()){
                        genegrouplist1.put(genei, genelist1);
                    }
                    
                    for(int xx = 0; xx < genelist0.size(); xx++){
                        List<Integer> genelist3 = new ArrayList<>();
                        int numi = genelist0.get(xx);
                        
                        if(genelist2.contains(numi))continue;
                        for(int yy = xx+1; yy < genelist0.size(); yy++){
                            int numj = genelist0.get(yy);
                            if(gene_cluster[numi][ew] == gene_cluster[numj][ew]){
                                genelist3.add(numj);
                                genelist2.add(numj);
                            }
                        }
                        if(!genelist3.isEmpty()){
                            genegrouplist1.put(numi, genelist3);
                        }
                    }
                }
                if(genegrouplist1.isEmpty())break;
                genemap.put(sw+"_"+ew, genegrouplist1);
            }           
        }
        System.out.println("group and combination succeed");
        
        //list all the gene_group
        int geneid = 0;
        for(int sw = 0; sw < endwindow - startwindow; sw++){
            for(int ew = sw + 1; ew <= endwindow - startwindow; ew++){
                HashMap<Integer, List<Integer>>genegrouplist = genemap.get(sw + "_" + ew);
                if(genegrouplist == null || genegrouplist.isEmpty())break;
                for (Map.Entry<Integer, List<Integer>> entry : genegrouplist.entrySet()) {
                    gene_group genegroup0 = new gene_group();
                    int gene0 = entry.getKey();
                    List<Integer> genelist = entry.getValue();
                    if(genelist.isEmpty())continue;
                    genegroup0.genelist.add(gene0);
                    genegroup0.genelist.addAll(genelist);
                    genegroup0.startw = sw;
                    genegroup0.endw = ew;
                    genegroupIDold.put(geneid, genegroup0);
                    geneid++;
                    
                    
                    fosg.write((geneid + ":\t").getBytes());
                    List<Integer> genelists = genegroup0.genelist;
                    for (Iterator<Integer> it0 = genelists.iterator(); it0.hasNext();) {
                        int genenum0 = it0.next();
                        fosg.write((genenum0 + " ").getBytes());
                    }
                     fosg.write(("\n\tstartwindow:\t" + genegroup0.startw + "\tendwindow\t" + genegroup0.endw + "\n").getBytes());
                }              
            }
        }
        System.out.println("list succeed\n" + geneid);
        
        //delete unnecessary and write table
        int geneid0 = 0;
        for(int i = 0; i < geneid; i++){
  
            boolean flag = true;
            gene_group genei = genegroupIDold.get(i);
            if(genei.genelist.size() < 4)continue;
            for(int j = 0; j < geneid; j++){
                if(i == j)continue;
                gene_group genej = genegroupIDold.get(j);
                if(genej.startw > genei.startw)break;
                if(genej.endw < genei.endw)continue;
                List<Integer> genelist1 = genei.genelist;
                List<Integer> genelist2 = genej.genelist;
                if(genelist2.containsAll(genelist1) && genelist1.size() == genelist2.size()){
                    flag = false;break;
                }             
            }
            if(flag){
                genegroupID.put(geneid0, genei);
                genegroup_trees.addVertex(geneid0);
                fost.write((geneid0 + ":\t").getBytes());
                geneid0++;
                List<Integer> genelists = genei.genelist;
                for (Iterator<Integer> it0 = genelists.iterator(); it0.hasNext();) {
                    int genenum0 = it0.next();
                    fost.write((genenum0 + " ").getBytes());
                }
                fost.write(("\n\tstartwindow:\t" + genei.startw + "\tendwindow\t" + genei.endw + "\n").getBytes());
            }
        }
        System.out.println("delete and write table succeed");
        System.out.println("id: " + geneid + "\tid0: \t" + geneid0);
        
        //connect nodes
        int edgenum = 0;
        int genegnumdiffw[] = new int[endwindow - startwindow + 2];
        int genegnumdiffgn[] = new int[200+1];
        for(int i = 2; i <= endwindow - startwindow+1; i++)genegnumdiffw[i] = 0;
        for(int i = 2; i <= 200; i++)genegnumdiffgn[i] = 0;
        for (int i = 0; i < geneid0; i++) {
            gene_group genei = genegroupID.get(i);
            fosxy.write((genei.genelist.size() + "," + (genei.endw - genei.startw + 1) + "\n").getBytes());
            genegnumdiffw[genei.endw - genei.startw+1]++;
            genegnumdiffgn[genei.genelist.size()]++;
            for (int j = 0; j < geneid0; j++) {
                if (i == j) {
                    continue;
                }
                gene_group genej = genegroupID.get(j);
                if (genej.startw > genei.startw) {
                    break;
                }
                if (genej.endw < genei.endw) {
                    continue;
                }
                List<Integer> genelist1 = genei.genelist;
                List<Integer> genelist2 = genej.genelist;
                if (genelist1.containsAll(genelist2)) {
                    List<Integer> childrenlist = children.get(i);
                    fosc.write((i + "\t" + j + "\n").getBytes());
                    //genegroup_trees.addEdge(i + "\t" + j, i, j);
                    if(childrenlist == null){
                        childrenlist = new ArrayList<>();
                    }
                    childrenlist.add(j);
                    for(Iterator<Integer> its = childrenlist.iterator(); its.hasNext();){
                        int c = its.next();
                        List<Integer> parentlist = parent.get(c);
                        if(parentlist == null)parentlist = new ArrayList<>();
                        parentlist.add(i);
                        parent.put(c, parentlist);
                    }
                    children.put(i, childrenlist);
                    edgenum++;
                }
            }
        }
        
        //delete unnecessary edge;
        int edgenum0 = 0;
        for(int i = 0; i < geneid0; i++){
            List<Integer> childrenlist = children.get(i);
            if(childrenlist == null || childrenlist.isEmpty()){
                fosnode.write((i + "\t" + 2 + "\n").getBytes());continue;
            }
            if(parent.get(i) == null || parent.get(i).isEmpty()){
                fosnode.write((i + "\t" + 1 + "\n").getBytes());
            }
            else{
                fosnode.write((i + "\t" + 0 + "\n").getBytes());
            }
            for(Iterator<Integer> it = childrenlist.iterator(); it.hasNext();){
                int gene0 = it.next();
                boolean flag = true;
                for(Iterator<Integer> it2 = childrenlist.iterator(); it2.hasNext();){
                    int gene1 = it2.next();
                    if(gene0 == gene1)continue;
                    if(children.get(gene1) == null)continue;
                    if(children.get(gene1).contains(gene0)){
                        flag = false;break;
                    }
                }
                if(flag){
                    foscn.write((i + "\tof\t" + gene0 + "\n").getBytes());
                    genegroup_trees.addEdge(i + "\t" + gene0, i, gene0);
                    edgenum0++;
                }
            }
        }
        System.out.println("edgenum: " + edgenum + "\tedgenum0: " + edgenum0);
        int s1 = 0, s2 = 0;
        for(int i = 2; i <= endwindow - startwindow+1; i++){
            System.out.println(genegnumdiffw[i]);
            s1 += genegnumdiffw[i];
        }
        System.out.println(s1 + " \t" + s2);
        for(int i = 2; i <= genenum; i++){
            System.out.println(genegnumdiffgn[i]);
            s2 += genegnumdiffgn[i];
        }
        
    }
}
