import org.apache.commons.math3.linear.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Scanner;
import javax.imageio.ImageIO;

public class SVD_Compression {
//*****************
	
	//converting image pixel RGB into grayscale
	public static double[][] grayscale(int height, int weight, BufferedImage img) {
		double[][] pixels = new double[height][weight];
		for(int i=0; i<height; i++)
            for (int j=0; j<weight; j++) {
                int rgb = img.getRGB(j, i),
                	r = (rgb >> 16) & 0xff, 
                	g = (rgb >> 8) & 0xff, 
                	b = (rgb & 0xff);
                pixels[i][j] = (r+g+b)/3;
            }
		return pixels;
	}
	
	//spits out Sigma by taking root of eigenvalues of AAt (and AtA) and (by reference) number of positive values
	public static RealMatrix getS(double[] eigenvalueMatrix, int height, int weight, int[] positive) {
		//creating S matrix by taking root of eigenvaluesMatrix
        double[][] S = new double[height][weight];
        positive[0]=0;
        for(int i=0; i<height; i++) {
        	for(int j=0; j<i; j++)
        		S[i][j] = 0;
        	
        	if(eigenvalueMatrix[i]>0) {
        		S[i][i] = Math.sqrt(eigenvalueMatrix[i]);
        		positive[0]++;
        	} else
        		S[i][i] = 0;
        	
        	for(int j=i+1; j<weight; j++)
        		S[i][j] = 0;
        }
        return MatrixUtils.createRealMatrix(S);
	}
		
//*****************
	
	public static void main(String[] args) {
//***************
		
	//*begin* *input elements*
		Scanner reader = new Scanner(System.in);
		
			//get image
        BufferedImage img = null;
        try {
    		System.out.println("insert path to image (only type name of file to point to root). for example:\nC:/Users/Hossein/Desktop/input.jpg");
        	String fileLocation = reader.next();
            img = ImageIO.read(new File(fileLocation));
        } catch (IOException e) {
            System.out.println("not found");
        }
        
        	//getting width and height
        int n = img.getWidth(), m = img.getHeight();
		
        	//converting image pixel RGB into grayscale
        double[][] pixels = grayscale(m,n,img);
        
			//getting compression rate and number of final eigenvalues
		System.out.println("compression ratio (any number between 0-100):");
        double perc_reduce = reader.nextInt();
        
	//*end* *input elements*

	//*begin* *SVD*
        //creating matrix of AAt and AtA
		RealMatrix image = MatrixUtils.createRealMatrix(pixels), imageT = image.transpose(), AAt = image.multiply(imageT), AtA = imageT.multiply(image);

		
		//eigenvalues and eigenDecomposition of AAt and AtA
	    EigenDecomposition eigenDeco_AAt = new EigenDecomposition(AAt), eigenDeco_AtA = new EigenDecomposition(AtA);
	    double[] eigenvalueMatrix = eigenDeco_AAt.getRealEigenvalues();
        
        //creating Sigma and getting number of positive values
	    int[] positive = new int[1];
        RealMatrix S = getS(eigenvalueMatrix, m, n, positive);
        
        //creating U
        RealMatrix U = MatrixUtils.createRealMatrix(m,m);
        for(int i=0; i<m; i++)
            U.setColumn(i,eigenDeco_AAt.getEigenvector(i).unitVector().toArray());

        //creating Vt
        RealMatrix Vt = MatrixUtils.createRealMatrix(n,n);
        for(int i=0; i<n; i++)
            Vt.setRow(i,eigenDeco_AtA.getEigenvector(i).unitVector().toArray());
    //*end* *SVD*

    //*begin* *creating compressed U,S,Vt*
        //k is percent of nonzero eigenvalues
        int k=(int) ((100-perc_reduce)*positive[0]*n/(100*(m+n+1)));
        RealMatrix newS = S.subtract(S);
        for(int i=0; i<k; i++)
            newS.setColumnMatrix(i,S.getColumnMatrix(i));
        
/*
        RealMatrix newU = U.subtract(U);
        for(int i=0; i<k; i++)
            newU.setRowMatrix(i,U.getRowMatrix(i));
        RealMatrix newVt = Vt.subtract(Vt);
        for(int i=0; i<k; i++)
            newVt.setColumnMatrix(i,Vt.getColumnMatrix(i));
*/
//
        RealMatrix newU = U.multiply(newS);
        RealMatrix newVt = newS.multiply(Vt);
        
        for(int i=0; i<k; i++)
            newS.setEntry(i,i,1/newS.getEntry(i,i));
        newS=newS.transpose();
//        
        RealMatrix newImage = newU.multiply(newS.multiply(newVt));
        
    //*end* *creating compressed U,S,Vt*
        
    //*begin* *Output*
        //creating new Buffer image
        BufferedImage newImg = new BufferedImage(n, m, BufferedImage.TYPE_BYTE_GRAY);
        for (int i=0; i<m; i++)
            for (int j=0; j < n; j++) {
                int value = (int)newImage.getEntry(i,j) << 16 | (int)newImage.getEntry(i,j) << 8 | (int)newImage.getEntry(i,j);
                newImg.setRGB(j, i, value);
            }
        
        //output image file
		System.out.println("insert path to output (only type name of file to output to root). for example:\nC:/Users/Hossein/Desktop/output.jpg");
    	String fileLocation = reader.next();
        File outputfile = new File(fileLocation);
        try {
            ImageIO.write(newImg, "jpg", outputfile);
        } catch (IOException e1) {}
        
        //output newU, newVt, newS
        int l=fileLocation.length()-1;
        while(l>=0 & fileLocation.charAt(l)!='/')
        	l--;
        fileLocation=fileLocation.substring(0, l+1);
        
        String text = "[";
        for(int i=0;i<k;i++)
        	text+=", "+S.getEntry(i, i);
        try (PrintWriter out = new PrintWriter(fileLocation+"newS.txt")) {
            out.println(text+"]");
        } catch (IOException e1) {}

        text = "";
        for(int i=0; i<k; i++)
        	text+=Arrays.toString(newVt.getRow(i))+" :nextline: ";
        try (PrintWriter out = new PrintWriter(fileLocation+"newVt.txt")) {
            out.println(text);
        } catch (IOException e1) {}
        
        text="";
        if(k==0) k++;
        for(int i=0; i<m; i++)
        	text+=Arrays.toString(Arrays.copyOfRange(newU.getRow(i), 0, k-1))+" :nextline: ";
        try (PrintWriter out = new PrintWriter(fileLocation+"newU.txt")) {
            out.println(text);
        } catch (IOException e1) {}

        //finishing message
        System.out.println("\ndone");
        
		reader.close();
     //*end* *Output*
        
//***************
	}
}
