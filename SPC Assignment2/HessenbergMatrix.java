import java.util.Random;

public class HessenbergMatrix {

	private static double matrix[][] = new double[1000][1000];
	private static double multiplier[][] = new double[1000][1000];
	private static double result[] = new double[1000];
	private static HessenbergMatrix hs = new HessenbergMatrix();
	
	public static void main(String[] args) {
		
		long startTime, stopTime, totalTime;
		startTime = System.currentTimeMillis();
		
		hs.matrixGenerator();
		hs.decomposition();
		hs.forwardSolving();
		hs.backwardSolving();
		
		stopTime = System.currentTimeMillis();
	    totalTime = stopTime - startTime;
	    System.out.println("Total Time Taken: " + totalTime + " ms");
	}
	
	private void matrixGenerator() {	
		Random rand = new Random();
		for (int i = 0; i < 1000; i++) {
			Integer constant = rand.nextInt()% 100;
			result[i] = Math.abs(constant);
			
	        for (int j = 0; j < 1000; j++) {
	        	if(i > (j + 1)){
	        		matrix[i][j] = 0;
	        	} else{
	        		Integer num = rand.nextInt()% 100; 
	        		matrix[i][j] = Math.abs(num);
	        	}
	        }
	    }	
	}
	
	private void decomposition() {
		for(int k = 0; k < 999; k++){
			for(int i = k+1; i < 1000; i++){
				multiplier[i][k] = matrix[i][k] / matrix[k][k];
				for(int j = k; j < 1000; j++){
					matrix[i][j] = matrix[i][j] - (multiplier[i][k] * matrix[k][j]);
				}
			}
		}	
		for(int i = 0; i < 1000; i++) {
			for(int j = 0; j < 1000; j++) {
				if(i == j){
					multiplier[i][j] = 1;
				}
			}
		}
	}
	
	private void forwardSolving() {	
		result[0] = result[0]/multiplier[0][0];
		for(int k = 1;k <= 999; k++){
			for(int j = k-1; j >= 0; j--){
				result[k] -= multiplier[k][j] * result[j];
			}
			result[k] = result[k]/multiplier[k][k];
		}
	}
	
	private void backwardSolving() {
		 result[999] = result[999]/matrix[999][999];
		 for(int k = 998; k >= 0; k--) {
		    for(int j = k+1; j < 1000; j++){
		    	result[k] -= matrix[k][j]*result[j];
		    }
		    result[k] = result[k]/matrix[k][k];
		 }
	}
}
