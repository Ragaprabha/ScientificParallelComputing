public class HessenbergResult {
	private static double matrix[][] = {{3,4,1,7}, {1,4,2,3}, {0,2,3,4}, {0,0,1,3}};
	private static double multiplier[][] = new double[4][4];
	private static double result[] = {1,2,3,4};
	private static HessenbergResult hs = new HessenbergResult();
	
	public static void main(String[] args) {
		
		long startTime, stopTime, totalTime;
		startTime = System.currentTimeMillis();
		
		System.out.print("The Upper Hessenberg Matrix is:\n");		
        for(int i=0; i<4; i++) {
			for(int j=0; j<4; j++) {
				System.out.print("\t" + matrix[i][j]);
			}
			System.out.println();
		}
        
        System.out.print("The constant vector is:\n");
        for(int i=0; i<4; i++) {
        	System.out.print("\t" + result[i]);
			System.out.println();
		}
        
		hs.decomposition();
		System.out.println("\nThe Upper Triangular Matrix is:");
		for(int i=0; i<4; i++) {
			for(int j=0; j<4; j++) {
				System.out.print("\t" + matrix[i][j]);
			}
			System.out.println();
		}
		System.out.println("\nThe Lower Triangular Matrix is:");
		for(int i=0; i<4; i++) {
			for(int j=0; j<4; j++) {
				System.out.print("\t" + multiplier[i][j]);
			}
			System.out.println();
		}
		
		hs.forwardSolving();
		System.out.println("\nThe Y Vector is:");
		
		for(int i=0; i<4; i++) {
			
				System.out.print("\t" + result[i]);
				System.out.println();
			}
		
		hs.backwardSolving();
		System.out.println("\nThe X Vector is:");
		
		for(int i=0; i<4; i++) {
			
				System.out.print("\t" + result[i]);
				System.out.println();
			}
		
		stopTime = System.currentTimeMillis();
	    totalTime = stopTime - startTime;
	    System.out.println("Total Time Taken: " + totalTime + " ms");
	}
	
	private void decomposition() {
		for(int k = 0; k < 3; k++){
			for(int i = k+1; i < 4; i++){
				multiplier[i][k] = matrix[i][k] / matrix[k][k];
				for(int j = k; j < 4; j++){
					matrix[i][j] = matrix[i][j] - (multiplier[i][k] * matrix[k][j]);
				}
			}
		}	
		for(int i = 0; i < 4; i++) {
			for(int j = 0; j < 4; j++) {
				if(i == j){
					multiplier[i][j] = 1;
				}
			}
		}
	}
	
	private void forwardSolving() {	
		result[0] = result[0]/multiplier[0][0];
		for(int k = 1;k <= 3; k++){
			for(int j = k-1; j >= 0; j--){
				result[k] -= multiplier[k][j] * result[j];
			}
			result[k] = result[k]/multiplier[k][k];
		}
	}
	
	private void backwardSolving() {
		 result[3] = result[3]/matrix[3][3];
		 for(int k = 2; k >= 0; k--) {
		    for(int j = k+1; j < 4; j++){
		    	result[k] -= matrix[k][j]*result[j];
		    }
		    result[k] = result[k]/matrix[k][k];
		 }
	}


}
