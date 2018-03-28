package stms;

import gnu.trove.map.hash.TIntIntHashMap;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map.Entry;

import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.random.JDKRandomGenerator;
import org.apache.commons.math3.special.Gamma;

import util.IOUtils;

public class SparsePTM {
	
	public int K1 = 500;
	public int K2 = 20;
	
	public int M;
	public int V;
	
	public double alpha = 0.1;
	
	public double alpha0 = 1E-12;
	public double alpha1 = 0.1;
	
	public double beta = 0.01;
	
	public double gamma0 = 1.0;
	public double gamma1 = 1.0;
	
	public boolean a[][];
	public int a_sum[];
	
	public double pi_a[];
	
	public JDKRandomGenerator rand;
	
	public int mp[];
	
	public int npk[][];
	public int npkSum[];
	
	public int nkw[][];
	public int nkwSum[];
	
	public int zAssigns_1[];
	public int zAssigns_2[][];
	
	public int niters = 500;
	
	public List<List<Integer>> docs = new ArrayList<List<Integer>>();
	public HashMap<String, Integer> w2i = new HashMap<String, Integer>();
	public HashMap<Integer, String> i2w = new HashMap<Integer, String>();
	
	public void loadTxts(String txtPath) {
		BufferedReader reader = IOUtils.getReader(txtPath, "UTF-8");
		
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
	
	public void initSTPAM() {
		rand = new JDKRandomGenerator();
		rand.setSeed(System.currentTimeMillis());
		
		BetaDistribution betaDist = new BetaDistribution(rand, gamma1 , gamma0);
		
		mp = new int[K1];
		
		npk = new int[K1][K2];
		npkSum = new int[K1];
		
		nkw = new int[K2][V];
		nkwSum = new int[K2];
		
		a = new boolean[K1][K2];
		a_sum = new int[K1];
		
		pi_a = new double[K1];
		
		for (int k1 = 0; k1 != K1; k1++) {
			pi_a[k1] = betaDist.sample();
			for (int k2 = 0; k2 != K2; k2++) {
				a[k1][k2] = true;
			}
			a_sum[k1] = K2;
		}
		
		zAssigns_1 = new int[M];
		zAssigns_2 = new int[M][];
		
		for (int m = 0; m != M; m++) {
			int N = docs.get(m).size();
			zAssigns_2[m] = new int[N];
			
			int z1 = (int) Math.floor(rand.nextDouble()*K1);
			zAssigns_1[m] = z1;
			
			mp[z1] ++;
			
			for (int n = 0; n != N; n++) {
				int w = docs.get(m).get(n);
				int z2 = (int) Math.floor(rand.nextDouble()*K2);
				
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
		
		TIntIntHashMap k2Count = new TIntIntHashMap();
		for (int n = 0; n != N; n++){
			int z2 = zAssigns_2[m][n];
			k2Count.adjustOrPutValue(z2, 1, 1);
			
			npk[z1][z2] --;
			npkSum[z1] --;
		}
		
		
		double[] pTable = new double[K1];
		
		for (int k = 0; k != K1; k++) {
			double expectTM = 1.0;
			int index = 0;
			for (int z2 : k2Count.keys()) {
				int x = a[k][z2] ? 1 : 0;
				int c = k2Count.get(z2);
				for (int i = 0; i != c; i++) {
					expectTM *= (npk[k][z2] + x*alpha1 + alpha0 + i) / (K2*alpha0 + a_sum[k]*alpha1 + npkSum[k] + index);
					index ++;
				}
			}
			
			pTable[k] = (mp[k] + alpha) / (M + K1 * alpha) * expectTM;
		}
		
		for (int k = 1; k != K1; k++) {
			pTable[k] += pTable[k-1];
		}
		
		double r = rand.nextDouble() * pTable[K1-1];
		
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
		
		double[] pTable = new double[K2];
		
		for (int k = 0; k != K2; k++) {
			int x = a[z1][k] ? 1 : 0;
			pTable[k] = (npk[z1][k] + x*alpha1 + alpha0) / (npkSum[z1] + K2*alpha0 + a_sum[z1]*alpha1) *
					(nkw[k][w] + beta) / (nkwSum[k] + VBeta);
		}
		
		for (int k = 1; k != K2; k++) {
			pTable[k] += pTable[k-1];
		}
		
		double r = rand.nextDouble() * pTable[K2-1];
		
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
	
	public void sampleBinaryAMatrix() {
		int GIBBS_ITER = 1;
		
		a_sum = new int[K1];
		for (int m = 0; m != K1; m++) {
			for (int k = 0; k != K2; k++) {
				a[m][k] = npk[m][k] > 0;
				a_sum[m] += a[m][k] ? 1 : 0;
			}
		}
		
		double log_diff, ratio, p;
		for (int iter = 0; iter != GIBBS_ITER; iter++) {
			for (int m = 0; m != K1; m++) {
				for (int k = 0; k != K2; k++) {
					if (a[m][k] && npk[m][k] == 0) {
						log_diff = Gamma.logGamma(a_sum[m]*alpha1 + K2*alpha0)
								- Gamma.logGamma((a_sum[m]-1)*alpha1 + K2*alpha0);
						log_diff -= Gamma.logGamma(npkSum[m] + a_sum[m]*alpha1 + K2*alpha0)
								- Gamma.logGamma(npkSum[m] + (a_sum[m]-1)*alpha1 + K2*alpha0);
						
						ratio = Math.exp(log_diff) * pi_a[m] / (1.0-pi_a[m]);
						p = ratio / (1.0 + ratio);
						if (rand.nextDouble() > p) { 
							a[m][k] = false;
							a_sum[m] --;
						}
					} else if (!a[m][k]) {
						log_diff = Gamma.logGamma((a_sum[m]+1)*alpha1 + K2*alpha0)
								- Gamma.logGamma(a_sum[m]*alpha1 + K2*alpha0);
						log_diff -= Gamma.logGamma(npkSum[m] + (a_sum[m]+1)*alpha1 + K2*alpha0)
								- Gamma.logGamma(npkSum[m] + a_sum[m]*alpha1 + K2*alpha0);
						
						ratio = Math.exp(log_diff) * pi_a[m] / (1.0-pi_a[m]);
						p = ratio / (1.0 + ratio);
						if (rand.nextDouble() < p) { 
							a[m][k] = true;
							a_sum[m] ++;
						}
					}
				}
				
				BetaDistribution betaDist = new BetaDistribution(rand, gamma1 + a_sum[m], gamma0 + K2 - a_sum[m]);
				pi_a[m] = betaDist.sample();
			}
		}
	}
	
	public void estimate() {
		for (int iter = 0; iter != niters; iter++) {
			System.out.println("SPTM Iteration: " + iter + " ...");
			long begin = System.currentTimeMillis();
			for (int m = 0; m != M; m++) {
				this.sampleZ1(m);
			}
			for (int m = 0; m != M; m++) {
				int N = docs.get(m).size();
				for (int n = 0; n != N; n++) {
					sampleZ2(m, n);
				}
			}
			
			if(iter % 1 == 0) {
				sampleBinaryAMatrix();
				System.out.print("\nPdoc sparsity : ");
				double sparsity = 0.0;
				for (int m = 0; m != K1; m++) {
					sparsity += 1.0 - pi_a[m];
				}
				System.out.println((int)(100*sparsity/K1));
				
				System.out.print("Pdoc real sparsity : ");
				sparsity = 0.0;
				for (int m = 0; m != K1; m++) {
					sparsity += 1.0 - (double)a_sum[m]/K2;
				}
				System.out.println((int)(100*sparsity/K1)+"\n");
			}
			System.out.println("Time(s) : "+(System.currentTimeMillis() - begin) / 1000.0+"\n");
		}
		return;	
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
				theta[m][k] = (ndk[m][k] + alpha) / (ndkSum[m] + K2*alpha);
			}
		}
		return theta;
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

	public void printTopics(int top_n) {
		double[][] phi = computePhi();
		ArrayList<List<Entry<String, Double>>> pairsList = this
				.sortedTopicWords(phi, K2);
		for (int k = 0; k != K2; k++) {
			System.out.println("Topic " + k + ":");
			for (int i = 0; i != top_n; i++) {
				System.out.println(pairsList.get(k).get(i).getKey() + " "
						+ pairsList.get(k).get(i).getValue());
			}
		}
	}
	
	public void savePhi(String path) {
		BufferedWriter writer = IOUtils.getWriter(path, "utf-8");
		
		double[][] phi = computePhi();
		
		try {
			for (int k = 0; k != K2; k++) {
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
	
	public void saveTheta(String path) {
		BufferedWriter writer = IOUtils.getWriter(path, "utf-8");
		
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
	
	public void saveWordmap(String path) {
		BufferedWriter writer = IOUtils.getWriter(path, "utf-8");
		
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
	
	public static void main(String args[]) {
		
		SparsePTM SPTM = new SparsePTM();
		SPTM.loadTxts("/Users/zuoyuan/Desktop/Experiment/20130601-origin");SPTM.initSTPAM();
		SPTM.estimate();
		SPTM.printTopics(10);
		
	}
}
