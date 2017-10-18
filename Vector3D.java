/** Class called Vector3D complete with properties representing the vector elements x,y,z.
* Getter and setter instance methods for the properties.
* A default constructor that sets all the elements to zero.
* A constructor vector that initializes the vector to (xx,yy,zz).
* A copy constructor for creating copies of vector object instances.
* Instance	methods that return	the	magnitude squared and the magnitude	of the vector.	
* Instance	method that	prints	out	the	elements of	the	vector	in	a useful format.	
* Instance	methods for	scalar	multiply and scalar	divide	by	a double.
* Static methods to	perform	vector addition, subtraction, cross product, and dot product operations on two Vectors.
* 
* @author JackFrankland s1404032
* @author Michael Spence s1419697
* @author Martin Tasker s1414001
* @version "10/2015"
*/

import java.util.Formatter;


public class Vector3D{

	//Create properties of vector, three components of 3D vector 
	
	private double xValue , yValue , zValue;

	/** Constructor Methods
	* Create three constructors 
	* Default constructor method to set all vector components to zero
	*/

	public Vector3D() {

		//Set to zero
		this.xValue = 0.0;
		this.yValue = 0.0;
		this.zValue = 0.0;
	}
	/** Constructor method to create vector given explicit components
	* @param xx a double giving the value of the x component of the vector
	* @param yy a double giving the value of the y component of the vector
	* @param zz a double giving the value of the z component of the vector
	*/

	public Vector3D(double xx, double yy, double zz) {

		this.xValue = xx;
		this.yValue = yy;
		this.zValue = zz;
	}

	/** Copy constructor to create create copies of Vector3D object instances
	* @param original the Vector3D to be copied
	*/

	public Vector3D(Vector3D original){

		this.setVector(original.getX() , original.getY() , original.getZ());



	}

	/** Instance Methods
	* Create getter and setter instance methods
	*/

	/** Create getter method
	* Getters provide access to private internal variables
  	*/
    

    /** Gets the x component  of a vector
    * @return double instance representing the x component of the vector
    */
	
	public double getX() {

		return this.xValue;
	} 

	/** Gets the y component  of a vector
    * @return double instance representing the y component of the vector
    */

	public double getY() {

		return this.yValue;
	}

	/** Gets the z component  of a vector
    * @return double instance representing the z component of the vector
    */

	public double getZ() {
		
		return this.zValue;

	}

	// Setters provide the internal variables

	/** Sets x component of vector
	* @param xx double to set the x component value
	*/

	public void setX(Double xx) {

		this.xValue = xx;
	}

	/** Sets y compoent of vector
	* @param yy double to set the y component value
	*/

	public void setY(Double yy) {

		this.yValue = yy;
	}

	/** Sets z compoent of vector
	* @param zz double to set the z component value
	*/

	public void setZ(Double zz) {

		this.zValue = zz;
	}

	/** Convenient setter method to set all three components 
	* @param xx a double to set x component value 
	* @param yy a double to set y component value 
	* @param zz a double to set z component value
	*/ 


	public void setVector(Double xx, Double yy, Double zz) {

		this.setX(xx);
		this.setY(yy);
		this.setZ(zz);
	}


	/** Instance method to return the magnitude squared of the vector
	* @return double instance representing the magnitude squared of the vector
	*/


	public double magVectorSquared() {

		return  this.getX()*this.getX() + this.getY()*this.getY() + this.getZ()*this.getZ();
	}

	/** Instance method to return the magnitude of the vector 
	* @return double instance representing the magnitude of the vector
	*/

	public double magVector() {

		return Math.sqrt(this.magVectorSquared());
	}

	/** Instance method tp print out components of vector in a useful string format
	* @return String representation of the vector instance
	*/ 

	public String toString() {

		double xx = this.getX();
		double yy = this.getY();
		double zz = this.getZ();

		Formatter f = new Formatter();

		return " ( " + f.format("%.3f , %.3f , %.3f" , xx , yy ,zz) + " ) ";
	}

	// Instance method for scalar multiplication and division by doubles

	/** Instance method for scalar multiplication 
	* @param a double scalar to multiply the vector by
	* @return Vector3D instance representing original vector multiplied by scalar
	*/

	public Vector3D scalarMult(Double a) {

		return new Vector3D(a*this.getX() , a*this.getY() , a*this.getZ());
	}

	/** Instance method for scalar division
	* @param a double scalar to divide the vector by
	* @return Vector3D instance representing original vector devided by a
	*/

	public Vector3D scalarDivison(Double a) {

		return new Vector3D( this.getX()/a , this.getY()/a , this.getZ()/a );
	}
	

	// Static Methods


	/** Here we use static methods so we can build compound statements 
	* such as a = Vector.addVector( b , Vector.subVector(c,d) );
	*/


	/** Static method for vector addition 
	* @param vectorA first vector to be added
	* @param vectorB second vector to be added
	* @return Vector3D instance representing vectorA + vectorB
	*/


	public static Vector3D addVector(Vector3D vectorA , Vector3D vectorB ) {

		return new Vector3D( vectorA.getX() + vectorB.getX() , vectorA.getY() + vectorB.getY() , vectorA.getZ() + vectorB.getZ() );


	}

	/** Static method for vector subtraction 
	* @param vectorA vector to be subtracted from
	* @param vectorB vector to be subtracted
	* @return Vector3D instance representing vectorA-vectorB
	*/


	public static Vector3D subVector(Vector3D vectorA , Vector3D vectorB ) {

		return new Vector3D( vectorA.getX() - vectorB.getX() , vectorA.getY() - vectorB.getY() , vectorA.getZ() - vectorB.getZ() );

	}

	/** Static method for calculating vector scalar product 
	 * @param vectorA first Vector to be in scalar prodcut
	 * @param vectorB second Vector to be in scalar product
	 * @return double instance represeting the scalar product of vectorA and vectorB
	 */

	public static double dotVector(Vector3D vectorA , Vector3D vectorB ) {

		return  vectorA.getX()*vectorB.getX() + vectorA.getY()*vectorB.getY() + vectorA.getZ()*vectorB.getZ();

	}

	/** Static method for calculating vector product
	* @param vectorA first vector in cross product
	* @param vectorB second vector in cross product
	* @return Vector3D instance representing the cross product of vectorA with vectorB
	*/

	public static Vector3D crossVector(Vector3D vectorA , Vector3D vectorB) {

		double xNew = (vectorA.getY() * vectorB.getZ()) - (vectorB.getY() * vectorA.getZ());
		double yNew = (vectorB.getX() * vectorA.getZ()) - (vectorA.getX() * vectorB.getZ());
		double zNew = (vectorA.getX() * vectorB.getY()) - (vectorB.getX() * vectorA.getY());

		return new Vector3D( xNew, yNew, zNew );
	}

	/** Static method vor determining of two vectors are the same
	* @param vectorA first vector being returned 
	* @param vectorB second vector being returned
	* @return Boolean true if they are the same, false if they are different
	*/

	public static Boolean equalVector(Vector3D vectorA, Vector3D vectorB) {

		if (vectorA==vectorB) {

			return true;
		}

		else {

			return false;
		}
	}





}