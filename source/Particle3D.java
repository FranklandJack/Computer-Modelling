	/** Class Called Particle3D.
* Vector	properties to	represent	the	position	and	velocity;
* Properties	to	represent	the	mass	and a text	label	for	the	particle;
* Getter	and	setter	methods	for	all the	properties
* A constructor	to	initialise	the	properties	with	specified	values;
* A default	constructor – choose	appropriate	defaults	for	the	properties;
* A toString method to	print	the	particle in a useful format
* A method to read	a	particle	from	a	Scanner object;
* A method	to	return	the	particle's	kinetic	energy;
* A	method	to	update	the	velocity of	the	particle	
* A	method	to	update	the	particle position	for	a	given timestep	and	the	current	particle velocity vector.
* A	method	to	update	the	particle position	for	a	given	timestep, current	particle	velocity	and	force.
* A	method	to	compute	the	relative separation	of	two	particles.
*
* @author Jack Frankland s1404032
* @author Michael Spence s1419697
* @author Martin Tasker  s1414001
* @version "23/02/2016"
*
**/

//We need the scanner class

import java.util.Scanner;

public class Particle3D {

/** Create properties of a 3D particle
* Position and velocity represented by vectors
* Mass represented by double, label represented by string
**/

private Vector3D position;
private Vector3D velocity;
private double mass;
private String label;

/** Constructor methods
* Create four Constructor methods
*/

/** Default constructor to set postition and velocity to zero vectors
* mass to zero double and label to a blank string 
*/

public Particle3D() {

	// Set to zero

	position = new Vector3D();
	velocity = new Vector3D();
	mass = 0.0;
	label = " ";

}

/** A second explicit constructor to initialize the properties with specific values
* @param p a Vector3D representing the position
* @param v a Vector3D representing the velocity 
* @param m a double representing the mass
* @param l a String representing the label
*/

public Particle3D(Vector3D p, Vector3D v, double m, String l) {

	this.setPosition(p);
	this.setVelocity(v);
	this.setMass(m);
	this.setLabel(l);
}

/** A third copy constructor to create copies of Particle3d object instnaces
* @param original a Particle3D instance to be copied
*/

public Particle3D(Particle3D original) {

	this.setParticle3D(original.getPosition(), original.getVelocity(), original.getMass(), original.getLabel());
}

/** Fourth constructor method that reads the values of the properties from a scanner object
* @param myScanner a Scanner from which the values of the properties are read
*/

public Particle3D(Scanner myScanner) {
	
	Vector3D position = new Vector3D(myScanner.nextDouble() , myScanner.nextDouble() , myScanner.nextDouble());

	Vector3D velocity = new Vector3D(myScanner.nextDouble() , myScanner.nextDouble() , myScanner.nextDouble());
	

	this.setPosition(position);
	this.setVelocity(velocity);
	this.setMass(myScanner.nextDouble());
	this.setLabel(myScanner.next());


}

// Setters and Getters

/** Get the position of the particle
* @return a Vector instance representing the position
*/

public Vector3D getPosition() {

	return position;
		
}

/** Get the velocity of the particle
* @return Vector instance represening the velocity 
*/

public Vector3D getVelocity() {

	return velocity; 
}

/** Get the mass of the particle 
* @return double instance representing the mass of the particle
*/

public double getMass() {

	return mass;
}

/** Get the label of the particle
* @return String instnace representing the label of the particle
*/

public String getLabel() {

	return label;
}


/** Set the position of a particle
*
* @param p a Vector represnting the position
*/

public void setPosition(Vector3D p) {

	this.position = p;

}

/** Set the velocity of a particle
*
* @param v a Vector representing the velocity 
*/

public void setVelocity(Vector3D v) {

	this.velocity = v;
}

/** Set the mass of a particle 
*
* @param m a double representing the mass
*/

public void setMass(double m) {

	this.mass = m;
}

/** Set the label of the particle
*
* @param l a Sting representing the label
*/

public void setLabel(String l) {

	this.label = l;

}

/** Convenient setter method to set all properties of Particle3D
* instance simulatnously
* @param p a Vector3D representing position
* @param v a Vector3D representing velocity
* @param m a double representing mass
* @param l a double representing label
*/

public void setParticle3D(Vector3D p, Vector3D v, double m, String l) {

	this.setPosition(p);
	this.setVelocity(v);
	this.setMass(m);
	this.setLabel(l);


}



//Instance Methods

/** toString method to print out properties of particle in useful String format
* @return String representation of the Particle3D instance
*/

public String toString() {

	double x = this.getPosition().getX();
	double y = this.getPosition().getY();
	double z = this.getPosition().getZ();
	String l  = this.getLabel();


	return l+ " " + x + " " +  y  + " " + z;
}


/** Instance method to calculate the kinetic energy of particle
* Kinetic energy calculated as 1/2*mv^2
* @return a double representing the kinetic energy of the particle
*/

public double kineticEnergy() {

	return 1.0/2 * this.getMass() * (this.getVelocity()).magVectorSquared();
}

/** A	method	to	update	the	velocity	of	the	particle	for	a	given timestep	and	force	vector,	as	
* v(t+dt) = v(t) + dt*f(t)/m
*
* @param dt a double that is the timestep
* @param force a Vector3D that is the current force on the particle
*/

public void leapVelocity(double dt, Vector3D force) {

	Vector3D timestepForceOverMass = force.scalarMult(dt/mass);

	velocity = Vector3D.addVector(velocity , timestepForceOverMass);

}

/** A	method	to	update	the	particle	position	for	a	given timestep	and	the	current	particle	
* velocity vector,	as	r(t+dt) = r(t) + dt*v(t);
*
* @param dt double that is the timestep
*/

public void leapPosition(double dt) {

	position = Vector3D.addVector(position , velocity.scalarMult(dt));

}

/** A	method	to	update	the	particle	position	for	a	given	timestep, current	particle	velocity	and	
* force,	as	r(t+dt) = r(t) + dt*v(t) +dt^2*f(t)/2m;
*
* @param dt a double that is the timestep
* @param force a Vector3D that is the current force on the particle
*/

public void leapPosition(double dt, Vector3D force) {

	Vector3D vectorA = force.scalarMult((dt*dt)/(2.0*mass));

	Vector3D vectorB = velocity.scalarMult(dt);

	Vector3D vectorC = Vector3D.addVector(vectorA , vectorB);



	position = Vector3D.addVector(position , vectorC);
}

/* An instnace method to create an array of forces between particles
* @param a Particle3D array instance from which the forces will be calculated
*/





//Static Methods

/* A static method to create an array of forces between particles
* @param a Particle3D array instance from which the forces will be calculated
* @param g double instance representing the gravational constant 
* @return a Vector3D array representing all the forces on the particles
*/

public static Vector3D[] forceArray(Particle3D[] a, double g) {
	int n = a.length;
	 Vector3D[] force = new Vector3D[n];
	
	for(int i = 0; i < n ; i++){
		Vector3D f = new Vector3D(0,0,0);
		
		for(int j = 0; j < n; j++){
			
			if(i != j){
			f = Vector3D.addVector(f, gravationalForce(a[i],a[j],g));
            }
				else{}
		}
         
        force[i] = f;
	}
	
	return force;
}

/** We use a static method to compute the relative seperation of particles
* so we can build complex statments.
* @param particleA a Particle3D representing the first particle 
* @param particleB a Particle3D representing the second particle
* @return deltaP a double representing the relative seperation of the particles
*/

public static Vector3D separation(Particle3D particleA, Particle3D particleB) {

	Vector3D separationVector = Vector3D.subVector(particleB.getPosition() , particleA.getPosition());


	return separationVector;
}

/** Static toString method to print out properties of particle in usefel String format
* This method acts on an array of particles taking the String representation of a specified element in the array
* @return String representation of Particle3D instance from an array
* @param the index of the element in the array you wish to print 
*/

public static String toString(int i, Particle3D[] a) {


	return a[i].toString();
}

/** Create methods to compute the gravational force as a vector and potential energy of the system as a double
* Based on the particle position and the mass at the origin
*/

/** Static method to compute the gravational force as a vector
* @param orbitingParticle a Particle3D instance representing the orbiting particle
* @param originParticle a Particle3D instance representing the mass at the origin
* @param g a double instnace representing the gravational constant 
* @return force a Vector3D instance representing the gravational force between the particles
*/

public static Vector3D gravationalForce(Particle3D orbitingParticle , Particle3D originParticle , double g) {

	double mass1 = orbitingParticle.getMass();
	double mass2 = originParticle.getMass();
	Vector3D separationVector = Particle3D.separation(originParticle,orbitingParticle);
	double separation = separationVector.magVector();
	
	
	return separationVector.scalarMult(-1.0*g* (mass1 * mass2)/Math.pow(separation , 3.0));


}


/** Static method to compute the potential energy of the system as a double 
* @param orbitingParticle a Particle3D instance representing the orbiting particle
* @param originParticle a Particle3D instance representing the mass at the origin
* @param g a double instance representing the gravational constant 
* @return potentialEnergy a double instance representing the potential energy of the system
*/

public static double gpEnergy(Particle3D orbitingParticle , Particle3D originParticle, double g) {

	double mass1 = orbitingParticle.getMass();
	double mass2 = originParticle.getMass();
	Vector3D separationVector = Particle3D.separation(originParticle,orbitingParticle);
	double separation = separationVector.magVector();
	
	
	return -1.0*g * mass1 * mass2 / separation;
}

/* a static leapPoisiton method that updates the position of a particles in an array
* @param dt  double instance time step for the velocity verlet algorithm
* @param a Particle3D array where the particles are held 
* @param f the array of vector forces acting on the particles
*/

public static void  leapPosition(double dt, Particle3D[] a, Vector3D[] f){
	
// for loop used here to update the position

for(int i = 0; i< a.length; i++){
	 a[i].leapPosition(dt, f[i]);
		}	
}

    /* a static leapPosition method that updates the position of particles in an array using the particle's velocity
     * @param dt double instance time step for the velocity verlet algorithm
     * @param a Particle3D array where particles are held
     */

 public static void leapPosition(double dt, Particle3D[] a) {
	
	for(int i = 0; i< a.length; i++){
	    a[i].leapPosition(dt);
	}
    }

/* a static method to calculate the total GP energy of an array of particles
* @param a the array of particles for which the GP energy will be calculated
* @param g double instance representing the gravational constant 
* @return double instance representing the total GP energy
*/

public static double totalPE(Particle3D[] a,double g){
	double p = 0;
	double n = a.length;
	for(int i=1; i<n; i++){
		for(int j = i+1; j<n; j++){
			p = p + gpEnergy(a[i], a[j],g);
		}

	}
	return p;
}

/* a static method to calculate the total kinetic energy of an array of particles
* @param a the array of particles for which the kinetic energy will be calculated
* @return double instance representing the total kinetic energy
*/
public static double totalKE(Particle3D[] a){
	double e = 0;
	double n = a.length;

	for(int i = 0; i<n; i++){
		e = e + a[i].kineticEnergy();

	}
	return e;
}


/* a static method that calculates the total energy of the system
* @param a Particle3D array which contains the system of particles
* @param g double instance representing the gravational constant 
* @return double instance representing the total energy of the system
*/

public static double totalEnergy(Particle3D[] a, double g) {
	return  Particle3D.totalPE(a,g) + Particle3D.totalKE(a);
}

/* a static method to update the velocities of the particles in the array
* @param System a Particle3D array instance in which the velocities will be updated 
* @param dt a double instance the time step for the algorithm 
* @param ForceArray a Vector3D array containing the forces acting on the particles
*/

public static void updateVelocityArray( Particle3D[] System, double dt, Vector3D[] ForceArray ) {
	    for (int i = 0; i < (System.length); i++) {
		System[i].leapVelocity(dt, ForceArray[i]);
	    }
    }

 /* a static method to read in the number of particles, and the properties of each particle
  * from a correctly formatted document
  */

public static Particle3D[] fillArray(Scanner particleScanner) {

	int particleNum = particleScanner.nextInt();

Particle3D[] pArray = new Particle3D[particleNum];

for(int i=0; i<pArray.length; i++) {

	pArray[i] = new Particle3D(particleScanner);

}
return pArray;
}
/* A static method to work out the angel between two particle on the x-y plane
 * with respect to the x-axis
 * @param p1 a Particle3D instance representing the first particle
 * @param p2 a Particle3D instance representing the second particle
 * @ return theta a double instance representing the angle
 */ 

public static double orbitAngle(Particle3D p1, Particle3D p2) {

	Vector3D d = Particle3D.separation(p1,p2);
	double y = d.getY();
	double x = d.getX();
	double theta = Math.atan(y/x);
	return theta;


}

}


