
public class Filter_Array_Computations {
	
	protected static int diam;
	protected static int[] BUFM;
	protected static int bufferWidth;
	protected static int SUMM;

	
	public static int[][] Filter2DArray (int[][] img2DArray, double radius, int filterType) {
		
		if (img2DArray == null) return null;
		
		int diam = (int) Math.ceil(radius*2);
		int filterSize = diam*diam;
		int width = img2DArray.length;
		int height = img2DArray[0].length;
		
		
		int[][] result = img2DArray;
		
		
		// initiate the box buffer
		int[] BUFM = new int[width];
		for (int i=0; i<width; i++) {
			for (int j=0; j<diam; j++) {
				BUFM[i] += img2DArray[i][j];
			}
		}
		
		float[] resultRow = new float[width];
		
		for (int h=diam; h<height-diam; h++) {
			if (h>diam) {
				for (int i=0; i<width; i++)
					BUFM[h] += (img2DArray[i][h+diam] - img2DArray[i][h]);
			}
			
			for (int i=0; i<diam; i++) {
				SUMM += BUFM[i];
			}
			
			result[diam][h] = (int)((float) SUMM / (float) filterSize);
			
			
			for (int w=diam+1; w<width-diam; w++) {
				
				SUMM += (BUFM[w+diam-1] - BUFM[w-1]);
				result[w][h] = (int)((float) SUMM / (float) filterSize);
				
			}	

		}
		
		return result;
		
	}
	
}
