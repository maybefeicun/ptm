package cn.edu.buaa.ptm;

/**
 * Created by dell on 2016/12/15.
 */
public class Main {
    public static void main(String args[]){
        boolean sparse = false;
        int P=0,K=0,iter=0;
        double alpha1=0,alpha2=0,beta=0;
        String doc_path="";
        if(args.length != 7&&args.length!=8){
            print_help();
            return;
        }
        int index = 0;
        if(args.length==8){
            if(!args[0].equals("sparse")){
                print_help();
                return;
            }
            index++;
            sparse = true;
        }
        boolean runs = true;
        try {
            P = Integer.parseInt(args[index++]);
            K = Integer.parseInt(args[index++]);
            iter = Integer.parseInt(args[index++]);
            alpha1 = Double.parseDouble(args[index++]);
            alpha2 = Double.parseDouble(args[index++]);
            beta = Double.parseDouble(args[index++]);
            doc_path = args[index++];
        }catch (Exception ex){
            print_help();
            runs = false;
        }

        if(runs){
            PseudoDocTM.PseudoDocTM(P, K, iter, 1, 1000, alpha1, alpha2, beta, 1, doc_path);
        }



    }

    public static void print_help(){
        System.out.println("Parameter error!");
        System.out.println("Usage: java -jar ptm.jar [sparse] P K iter alpha1 alpha2 beta doc_path");
        System.out.println("P\t:\tnumber of pesudo document.");
        System.out.println("K\t:\tnumber of topic.");
        System.out.println("iter\t:\tnumber of iteration times.");
        System.out.println("alpha1\t:\tprior parameter alpha1.");
        System.out.println("alpha2\t:\tprior parameter alpha2.");
        System.out.println("beta\t:\tprior parameter beta.");
        System.out.println("doc_path\t:\tpath of train file.(must be absolute path!)");
    }
}
