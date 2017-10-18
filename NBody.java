/** Main nethod for NBody simulation
* @author Jack Frankland s1404032
* @author Martin Tasker
* @author Michael Spence 
* @version 26/1/2016
*/
 
// Import relevent packages 
import java.util.Scanner; 
import java.io.*;
 
public class NBody{
     
    public static void main(String argv[]) throws IOException{
 
        double t = 0.0;
//
 
 
/** The main method will read in data from two files from the command line. 
* The fist argument in the command line will be the file containg the particle data.
* The second argument in the command line will be the file containing simulation perameters.
*/
 
// Create buffered reader to read in file text for the particle data.
 
BufferedReader particleFile = new BufferedReader(new FileReader(argv[0]));
 
// Attach a scanner to the file to parse in the particel data.
 
Scanner particleScanner = new Scanner(particleFile);
 
// Create buffered reader to read in file text for the parameters for the simulation.
 
BufferedReader simulationFile = new BufferedReader(new FileReader(argv[1]));
 
// Attach a scanner to the file to parse in the simulation parameters.
 
Scanner simulationScanner = new Scanner(simulationFile);
 

 
 
/* The particle data will all be stored in a 2D array
* We can use a method fromt the particle3D class to fill the array using data from the file
*/

 Particle3D[] pArray = Particle3D.fillArray(particleScanner);

// Read the simulation parameters into the file
double dt = simulationScanner.nextDouble();
int numsteps = simulationScanner.nextInt();
double gravConst = simulationScanner.nextDouble();
// an int to allow us to only print every tenth postition value to the file
int printchecker = 0;

//Loop to calculate total mass
double mass_tot = 0.0;
for (int q=0; q<pArray.length; q++) {
    double mass_p = pArray[q].getMass();
    mass_tot = mass_tot + mass_p;
}
//Loop to calculate total momentum
 Vector3D mom_tot = new Vector3D();
for (int r=0; r<pArray.length; r++) {
    Vector3D mom_p = pArray[r].getVelocity().scalarMult(pArray[r].getMass());
    mom_tot =  Vector3D.addVector(mom_tot, mom_p);
}
//Find compensation velocity
Vector3D v_comp = mom_tot.scalarDivison(mass_tot);

//Change initial velocity of each particle to initial v minus v_comp
for (int l=0; l<pArray.length; l++) {
    Vector3D newVelocity = Vector3D.subVector(pArray[l].getVelocity(), v_comp);
    pArray[l].setVelocity(newVelocity);
}

/* The main method will output the simulation results to a file named by the third argument in the command line
* Open the output file
* Attach a printwriter to it to output the results.
*/

String dataOutputFile = argv[2];
PrintWriter dataOutput = new PrintWriter(new FileWriter(dataOutputFile));

//Create the file to output the perhelion and aphelion data to and attach a print writer to it to output results
String periAphFile = "PerihelionAphelion.txt";
PrintWriter periAphOutput = new PrintWriter(new FileWriter(periAphFile));

//create the file to output the energy values and attach a print writer to it to output results
String energyFile = "Energy.txt";
PrintWriter energyOutput = new PrintWriter(new FileWriter(energyFile));

//Creat the file to output the orbit numbers of each particle
String orbitFile ="Orbits.txt";
PrintWriter orbitOutput = new PrintWriter(new FileWriter(orbitFile));

  // Create an array to hold the initial perihelion of the planets from the sun 
 double[] perihelionArray = new double [pArray.length];
for (int i = 1; i<pArray.length;i++) {
    Vector3D perihelionVector = Particle3D.separation(pArray[0],pArray[i]);
    perihelionArray[i] = perihelionVector.magVector();
}

  // Create an array to hold the initial aphelions of the planets from the sun 
 double[] aphelionArray = new double [pArray.length];
for (int i = 1; i<pArray.length;i++) {
    Vector3D aphelionVector = Particle3D.separation(pArray[0],pArray[i]);
    aphelionArray[i] = aphelionVector.magVector();
}

// Create an array to monitor orbits
int l = pArray.length;
double[] orbits = new double[l];
for(int i =1; i<l; i++){
orbits[i] = 0;
}

//Calculate the initial energy and output it to the energy file 

double totalEnergy = Particle3D.totalEnergy(pArray,gravConst);
energyOutput.printf("%g \n",totalEnergy);

 
/* We will create a loop to output the calculated values to a file
* First we will print the initial values.
*/
 
for(int i=0; i<numsteps; i++) {

    printchecker = printchecker + 1;
    int print = printchecker % 10;

if(print == 0){
//Print number of particles
    dataOutput.printf("%s \n" , pArray.length);
//Print current point value
    dataOutput.printf("point= %d \n" , i+1);
 }

 //Create array of old particle positions, before array is updated
    double[] oldY = new double[l];
    for(int g = 1; g<l; g++){
        Vector3D sep = Particle3D.separation(pArray[0], pArray[g]);
        oldY[g] = sep.getY();
    }



//update the perahelion for the 'ith' timestep for every particle

   for(int w=1;w<pArray.length;w++){
Vector3D x = Particle3D.separation(pArray[0], pArray[w]);
double xmag = x.magVector();
if(xmag<perihelionArray[w]){
    perihelionArray[w]=xmag;
}
else{}
}

//update the aphelion for the 'ith' timestep for every particle

   for(int u=1;u<pArray.length;u++){
Vector3D x = Particle3D.separation(pArray[0], pArray[u]);
double xmag = x.magVector();
if(xmag>aphelionArray[u]){
    aphelionArray[u]=xmag;
}
else{}
}

//Calculate the total energy at the 'ith' time step and output it to the energy file 

totalEnergy = Particle3D.totalEnergy(pArray,gravConst);
energyOutput.printf("%g \n",totalEnergy);

// the if statement means only every 10th value is printed to the file
    if(print==0){
//Nested loop to output the positions of the particles at point 'i'
    for (int j=0; j<pArray.length; j++){
        dataOutput.printf("%s \n" , Particle3D.toString(j, pArray));
        } 
    }
       
    Vector3D[] forceArray = new Vector3D[pArray.length];
    Vector3D[] forceArray_old = Particle3D.forceArray(pArray,gravConst);
    Particle3D.leapPosition(dt, pArray, forceArray_old);
    Vector3D[] forceArray_new = Particle3D.forceArray(pArray,gravConst);
    
    for (int k = 0; k < pArray.length; k++) {
            forceArray[k] = Vector3D.addVector(forceArray_old[k].scalarMult(0.5),forceArray_new[k].scalarMult(0.5));
    }
            Particle3D.updateVelocityArray(pArray,dt, forceArray);
            t = t + dt;

    /*We check to see if the x coordinate of every particle has changed signs. 
     *If it has then it tells it's crossed the x-z plane, and will be counted as half an orbit
     */

   for(int q =1;q<l;q++){
    Vector3D newSep = Particle3D.separation(pArray[0], pArray[q]);
    double newY = newSep.getY();
    if((oldY[q]*newY)<0){
        orbits[q]= orbits[q]+0.5;
    }
    else{
        orbits[q]=orbits[q];
    }
   }
  
}


//Print the perahelion and aphelion out to the file 

for(int w=1;w<pArray.length;w++){

     periAphOutput.printf("%s: aphelion %g , perihelion %g \n",pArray[w].getLabel() , aphelionArray[w] , perihelionArray[w]);
}

//Print the number of orbits for each particle in the simulation to a file

for(int w=1;w<pArray.length;w++){

    orbitOutput.printf("%s orbit: %g \n" , pArray[w].getLabel(), orbits[w]);
}
//close the output files
 periAphOutput.close();
 dataOutput.close();
 energyOutput.close();
 orbitOutput.close();    
}
 
 
}