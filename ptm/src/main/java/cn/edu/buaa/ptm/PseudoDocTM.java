package cn.edu.buaa.ptm;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

public class PseudoDocTM implements Runnable {

	public int K1 = 1000;
	public int K2 = 100;
	
	public int M;
	public int V;

	public double alpha1 = 0.1;
	public double alpha2 = 0.1;

	public double beta = 0.01;

	public int mp[];

	public int npk[][];
	public int npkSum[];

	public int nkw[][];
	public int nkwSum[];

	public int zAssigns_1[];
	public int zAssigns_2[][];

	public int niters = 200;
	public int saveStep = 1000;
	public String inputPath="";
	public String outputPath="";

	public int innerSteps = 10;

	public List<List<Integer>> docs = new ArrayList<List<Integer>>();
	public HashMap<String, Integer> w2i = new HashMap<String, Integer>();
	public HashMap<Integer, String> i2w = new HashMap<Integer, String>();


	public PseudoDocTM(int P,int K,int iter,int innerStep,int saveStep,double alpha1,double alpha2,double beta,String inputPath,String outputPath){
		this.K1=P;
		this.K2=K;
		this.niters=iter;
		this.innerSteps= innerStep;
		this.saveStep =saveStep;
		this.alpha1=alpha1;
		this.alpha2= alpha2;
		this.beta = beta;
		this.inputPath=inputPath;
		this.outputPath=outputPath;
	}

	public void loadTxts(String txtPath) {
		BufferedReader reader = IOUtil.getReader(txtPath, "UTF-8");

		String line;
		try {
			line = reader.readLine();
			while (line != null) {
				List<Integer> doc = new ArrayList<Integer>();

				String[] tokens = line.trim().split("\\s+");
				for (String token : tokens) {
					if (!w2i.containsKey(token)) {
						w2i.put(token, w2i.size());
						i2w.put(w2i.get(token), token);
					}
					doc.add(w2i.get(token));
				}
				docs.add(doc);
				line = reader.readLine();
			}
			reader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}


		M = docs.size();
		V = w2i.size();

		return;
	}

	public void initModel() {
		mp = new int[K1];

		npk = new int[K1][K2];
		npkSum = new int[K1];

		nkw = new int[K2][V];
		nkwSum = new int[K2];

		zAssigns_1 = new int[M];
		zAssigns_2 = new int[M][];

		for (int m = 0; m != M; m++) {
			int N = docs.get(m).size();
			zAssigns_2[m] = new int[N];

			int z1 = (int) Math.floor(Math.random()*K1);
			zAssigns_1[m] = z1;

			mp[z1] ++;

			for (int n = 0; n != N; n++) {
				int w = docs.get(m).get(n);
				int z2 = (int) Math.floor(Math.random()*K2);

				npk[z1][z2] ++;
				npkSum[z1] ++;

				nkw[z2][w] ++;
				nkwSum[z2] ++;

				zAssigns_2[m][n] = z2;
			}
		}
	}

	public void sampleZ1(int m) {
		int z1 = zAssigns_1[m];
		int N = docs.get(m).size();

		mp[z1] --;

		Map<Integer, Integer> k2Count = new HashMap<Integer, Integer>();
		for (int n = 0; n != N; n++){
			int z2 = zAssigns_2[m][n];
			if (k2Count.containsKey(z2)) {
				k2Count.put(z2, k2Count.get(z2)+1);
			} else {
				k2Count.put(z2, 1);
			}

			npk[z1][z2] --;
			npkSum[z1] --;
		}

		double k2Alpha2 = K2 * alpha2;

		double[] pTable = new double[K1];

		for (int k = 0; k != K1; k++) {
			double expectTM = 1.0;
			int index = 0;
			for (int z2 : k2Count.keySet()) {
				int c = k2Count.get(z2);
				for (int i = 0; i != c; i++) {
					expectTM *= (npk[k][z2] + alpha2 + i) / (k2Alpha2 + npkSum[k] + index);
					index ++;
				}
			}

			pTable[k] = (mp[k] + alpha1) / (M + K1 * alpha1) * expectTM;
		}

		for (int k = 1; k != K1; k++) {
			pTable[k] += pTable[k-1];
		}

		double r = Math.random() * pTable[K1-1];

		for (int k = 0; k != K1; k++) {
			if (pTable[k] > r) {
				z1 = k;
				break;
			}
		}

		mp[z1] ++;
		for (int n =0; n != N; n++) {
			int z2 = zAssigns_2[m][n];
			npk[z1][z2] ++;
			npkSum[z1] ++;
		}

		zAssigns_1[m] = z1;
	}
	
	public void sampleZ2(int m, int n) {
		int z1 = zAssigns_1[m];
		int z2 = zAssigns_2[m][n];
		int w = docs.get(m).get(n);

		npk[z1][z2] --;
		npkSum[z1] --;
		nkw[z2][w] --;
		nkwSum[z2] --;

		double VBeta = V * beta;
		double k2Alpha2 = K2 * alpha2;

		double[] pTable = new double[K2];

		for (int k = 0; k != K2; k++) {
			pTable[k] = (npk[z1][k] + alpha2) / (npkSum[z1] + k2Alpha2) *
					(nkw[k][w] + beta) / (nkwSum[k] + VBeta);
		}

		for (int k = 1; k != K2; k++) {
			pTable[k] += pTable[k-1];
		}

		double r = Math.random() * pTable[K2-1];

		for (int k = 0; k != K2; k++) {
			if (pTable[k] > r) {
				z2 = k;
				break;
			}
		}

		npk[z1][z2] ++;
		npkSum[z1] ++;
		nkw[z2][w] ++;
		nkwSum[z2] ++;

		zAssigns_2[m][n] = z2;
		return;
	}

	public void estimate() {
		long start = 0;
		for (int iter = 0; iter != niters; iter++) {
			start = System.currentTimeMillis();
			System.out.println("PAM4ST Iteration: " + iter + " ...");
			if(iter%this.saveStep==0&&iter!=0&&iter!=this.niters-1){
				this.storeResult(iter);
			}
			for (int i = 0; i != innerSteps; i++) {
				for (int m = 0; m != M; m++) {
					this.sampleZ1(m);
				}
			}
			for (int i = 0; i != innerSteps; i++) {
				for (int m = 0; m != M; m++) {
					int N = docs.get(m).size();
					for (int n = 0; n != N; n++) {
						sampleZ2(m, n);
					}
				}
			}
			System.out.println("cost time:"+(System.currentTimeMillis()-start));
		}
		return;
	}
	
	public double[][] computeThetaP() {
		double[][] theta = new double[K1][K2];
		for (int k1 = 0; k1 != K1; k1++) {
			for (int k2 = 0; k2 != K2; k2++) {
				theta[k1][k2] = (npk[k1][k2] + alpha2) / (npkSum[k1] + K2*alpha2);
			}
		}
		return theta;
	}
	
	public void saveThetaP(String path) throws IOException {
		BufferedWriter writer = IOUtil.getWriter(path);
		double[][] theta = this.computeThetaP();
		for (int k1 = 0; k1 != K1; k1++) {
			for (int k2 = 0; k2 != K2; k2++) {
				writer.append(theta[k1][k2]+" ");
			}
			writer.newLine();
		}
		writer.flush();
		writer.close();
	}
	
	public void saveZAssigns1(String path) throws IOException {
		BufferedWriter writer = IOUtil.getWriter(path);
		
		for (int m = 0; m != M; m++) {
			writer.append(zAssigns_1[m]+"\n");
		}
		
		writer.flush();
		writer.close();
	}

	public double[][] computePhi() {
		double[][] phi = new double[K2][V];
		for (int k = 0; k != K2; k++) {
			for (int v = 0; v != V; v++) {
				phi[k][v] = (nkw[k][v] + beta) / (nkwSum[k] + V*beta);
			}
		}
		return phi;
	}

	public ArrayList<List<Entry<String, Double>>> sortedTopicWords(
			double[][] phi, int T) {
		ArrayList<List<Entry<String, Double>>> res = new ArrayList<List<Entry<String, Double>>>();
		for (int k = 0; k != T; k++) {
			HashMap<String, Double> term2weight = new HashMap<String, Double>();
			for (String term : w2i.keySet())
				term2weight.put(term, phi[k][w2i.get(term)]);

			List<Entry<String, Double>> pairs = new ArrayList<Entry<String, Double>>(
					term2weight.entrySet());
			Collections.sort(pairs, new Comparator<Entry<String, Double>>() {
				public int compare(Entry<String, Double> o1,
						Entry<String, Double> o2) {
					return (o2.getValue().compareTo(o1.getValue()));
				}
			});
			res.add(pairs);
		}
		return res;
	}


	public void printTopics(String path,int top_n) throws IOException {
		BufferedWriter writer = IOUtil.getWriter(path);
		double[][] phi = computePhi();
		ArrayList<List<Entry<String, Double>>> pairsList = this
				.sortedTopicWords(phi, K2);
		for (int k = 0; k != K2; k++) {
			writer.write("Topic " + k + ":\n");
			for (int i = 0; i != top_n; i++) {
				writer.write(pairsList.get(k).get(i).getKey() + " "
						+ pairsList.get(k).get(i).getValue()+"\n");
			}
		}
		writer.close();
	}

	public void savePhi(String path) {
		BufferedWriter writer = IOUtil.getWriter(path, "utf-8");

		double[][] phi = computePhi();
		int K = phi.length;
		assert K > 0;
		int V = phi[0].length;

		try {
			for (int k = 0; k != K; k++) {
				for (int v = 0; v != V; v++) {
					writer.append(phi[k][v]+" ");
				}
				writer.append("\n");
			}
			writer.flush();
			writer.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return;
	}

	public void saveWordmap(String path) {
		BufferedWriter writer = IOUtil.getWriter(path, "utf-8");

		try {
			for (String word : w2i.keySet())
				writer.append(word + "\t" + w2i.get(word) + "\n");

			writer.flush();
			writer.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return;
	}

	public void saveAssign(String path){
		BufferedWriter writer = IOUtil.getWriter(path, "utf-8");
		try {
			for(int i=0;i<zAssigns_2.length;i++){
				for(int j=0;j<zAssigns_2[i].length;j++){
					writer.write(docs.get(i).get(j)+":"+zAssigns_2[i][j]+" ");
				}
				writer.write("\n");
			}
			writer.flush();
			writer.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return;
	}
	public void printModel(){
		System.out.println("\tK1 :"+this.K1+
				"\tK2 :"+this.K2+
				"\tniters :"+this.niters+
				"\tinnerSteps :"+this.innerSteps+
				"\tsaveStep :"+this.saveStep +
				"\talpha1 :"+this.alpha1+
				"\talpha2 :"+this.alpha2+
				"\tbeta :"+this.beta +
				"\tinputPath :"+this.inputPath+
				"\toutputPath :"+this.outputPath);
	}
	
	int[][] ndk;
	int[] ndkSum;
	
	public void convert_zassigns_to_arrays_theta(){
		ndk = new int[M][K2];
		ndkSum = new int[M];
		
		for (int m = 0; m != M; m++) {
			for (int n = 0; n != docs.get(m).size(); n++) {
				ndk[m][zAssigns_2[m][n]] ++;
				ndkSum[m] ++;
			}
		}
	}
	
	public double[][] computeTheta() {
		convert_zassigns_to_arrays_theta();
		double[][] theta = new double[M][K2];
		for (int m = 0; m != M; m++) {
			for (int k = 0; k != K2; k++) {
				theta[m][k] = (ndk[m][k] + alpha2) / (ndkSum[m] + K2 * alpha2);
			}
		}
		return theta;
	}
	
	public void saveTheta(String path) {
		BufferedWriter writer = IOUtil.getWriter(path, "utf-8");
		
		double[][] theta = computeTheta();
		try {
			for (int m = 0; m != M; m++) {
				for (int k = 0; k != K2; k++) {
					writer.append(theta[m][k]+" ");
				}
				writer.append("\n");
			}
			writer.flush();
			writer.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return;
	}
	
	public void storeResult(int times){
		String appendString="final";
		if(times!=0){
			appendString =times+"";
		}
		try {
			this.printTopics(outputPath+"/model-"+appendString+".twords",20);
			this.saveWordmap(outputPath+"/wordmap.txt");
			this.savePhi(outputPath+"/model-"+appendString+".phi");
			this.saveAssign(outputPath+"/model-"+appendString+".tassign");
			this.saveTheta(outputPath+"/model-"+appendString+".theta");
			this.saveThetaP(outputPath+"/model-"+appendString+".thetap");
			this.saveZAssigns1(outputPath+"/model-"+appendString+".assign1");
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	public void run() {
		printModel();
		this.loadTxts(inputPath);//
		this.initModel();
		this.estimate();
		this.storeResult(0);

	}
	
	
	public static void PseudoDocTM(int P,int K,int iter,int innerStep,int saveStep,double alpha1,double alpha2,double beta,int threadNum,String path){
		File trainFile = new File(path);
		String parent_path = trainFile.getParentFile().getAbsolutePath();
		(new File(parent_path+"/PTM_with_case_"+P+"_"+K+"_"+iter+"_"+alpha1+"_"+alpha2+"_"+beta+"/")).mkdirs();
		try {
			Thread.sleep(1000);
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
		(new PseudoDocTM(P,K,iter,innerStep,saveStep,alpha1,alpha2,beta,path,parent_path+"/PTM_with_case_"+P+"_"+K+"_"+iter+"_"+alpha1+"_"+alpha2+"_"+beta)).run();

	}
}
