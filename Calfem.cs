using System;
using System.Collections.Generic;
using System.Linq;
//using System.Text;
//using System.Collections;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

namespace Calfem
{
    public class ResultFEM
    {
        public Matrix<double> u;
        public Matrix<double> Q;

        public ResultFEM(Matrix<double> newU, Matrix<double> newQ)
        {
            u = newU;
            Q = newQ;
        }

        public void Print()
        {
            System.Console.WriteLine("Displacements");
            System.Console.WriteLine(u);
            System.Console.WriteLine("Boundary forces");
            System.Console.WriteLine(Q);
        }
    }

    public static class FEM
    {
        public static Matrix Spring1e(double ep)
        {
            //Element stiffnesmatrix for spring
            Matrix Ke = DenseMatrix.OfArray(new double[,] { { ep, -ep }, { -ep, ep } });
            return Ke;
        }

        public static Matrix<double> Assem(Matrix<double> K, Matrix<double> Ke, int[] edofRow)
        {
            for (int i = 0; i < edofRow.Length; i++)
            {
                for (int j = 0; j < edofRow.Length; j++)
                {
                    // Stiffness matrix assembly
                    K[edofRow[i], edofRow[j]] = K[edofRow[i], edofRow[j]] + Ke[i, j];
                }
            }
            return K;
        }

        public static void Assem(ref Matrix<double> K, Matrix<double> Ke, ref Matrix<double> F, Matrix<double> Fe, int[] edofRow)
        {
            for (int i = 0; i < edofRow.Length; i++)
            {
                for (int j = 0; j < edofRow.Length; j++)
                {
                    // Stiffness matrix assembly
                    K[edofRow[i], edofRow[j]] = K[edofRow[i], edofRow[j]] + Ke[i, j];
                }
                // Force matrix assembly
                F[edofRow[i], 0] = F[edofRow[i], 0] + Fe[i, 0];
            }
        }

        public static ResultFEM Solve(Matrix<double> K, Matrix<double> f, int[] bc, double[] bcVal)
        {
            // Create dof array:
            int ndof = f.RowCount;

            // Initialize displacement vector
            Matrix u = new DenseMatrix(ndof, 1);

            int[] dof = new int[ndof];
            for (int i = 0; i < ndof; i++)
            {
                dof[i] = i;
            }
            // Create the free dofs:
            int[] fdof = dof.Except(bc).ToArray();

            int nfdof = fdof.Length; // Number of free dofs.
            int nbdof = ndof - nfdof; //Constrained dofs

            // Create the permutation array which puts the constrained dofs last:
            int[] permute = bc.Union(fdof).ToArray();

            Permutation perm = new Permutation(permute);
            Permutation invPerm = perm.Inverse();
            Matrix<double> KPermuted = DenseMatrix.OfMatrix(K);
            KPermuted.PermuteRows(invPerm);
            KPermuted.PermuteColumns(invPerm);

            // Split K:::
            Matrix<double> Kff = KPermuted.SubMatrix(nbdof, nfdof, nbdof, nfdof);
            System.Console.WriteLine(Kff);  //Down right corner of matrix
            Matrix<double> Kfp = KPermuted.SubMatrix(nbdof, nfdof, 0, nbdof); //Down left corner of matrix

            // Split F:::
            Matrix<double> fPermuted = DenseMatrix.OfMatrix(f);
            fPermuted.PermuteRows(invPerm);
            Matrix<double> ff = fPermuted.SubMatrix(nbdof, nfdof, 0, 1);

            // Split displacements:::
            for (int i = 0; i < nbdof; i++)
            {
                u[i, 0] = bcVal[i];
            }

            Matrix<double> up = u.SubMatrix(0, nbdof, 0, 1); // Must set up to constrained values in bc.

            // Solve for the unknown displacements:::
            Matrix<double> s = Kff.Solve(ff.Subtract(Kfp.Multiply(up)));
            u.SetSubMatrix(nbdof, 0, s); // Set displacements back.

            System.Console.WriteLine(u);
            // Permute back u
            u.PermuteRows(perm);

            System.Console.WriteLine("U after permut");
            System.Console.WriteLine(K);
            System.Console.WriteLine(u);
            System.Console.WriteLine(f);

            // Get reaction forces:
            Matrix<double> Q = K.Multiply(u).Subtract(f);

            ResultFEM result = new ResultFEM(u, Q);

            return result;

        }

        public static ResultFEM Solve(Matrix<double> K, Matrix<double> f, int[] bc)
        {

            // Create dof array:
            int ndof = f.RowCount;

            // Initialize displacement vector
            Matrix u = new DenseMatrix(ndof, 1);

            int[] dof = new int[ndof];
            for (int i = 0; i < ndof; i++)
            {
                dof[i] = i;
            }
            // Create the free dofs:
            int[] fdof = dof.Except(bc).ToArray();

            int nfdof = fdof.Length; // Number of free dofs.
            int nbdof = ndof - nfdof; //Constrained dofs

            // Create the permutation array which puts the constrained dofs last:
            int[] permute = bc.Union(fdof).ToArray();

            Permutation perm = new Permutation(permute);
            Permutation invPerm = perm.Inverse();
            Matrix<double> KPermuted = DenseMatrix.OfMatrix(K);
            KPermuted.PermuteRows(invPerm);
            KPermuted.PermuteColumns(invPerm);

            // Split K:::
            Matrix<double> Kff = KPermuted.SubMatrix(nbdof, nfdof, nbdof, nfdof);
            System.Console.WriteLine(Kff);  //Down right corner of matrix
            Matrix<double> Kfp = KPermuted.SubMatrix(nbdof, nfdof, 0, nbdof); //Down left corner of matrix

            // Split F:::
            Matrix<double> fPermuted = DenseMatrix.OfMatrix(f);
            fPermuted.PermuteRows(invPerm);
            Matrix<double> ff = fPermuted.SubMatrix(nbdof, nfdof, 0, 1);


            Matrix<double> up = u.SubMatrix(0, nbdof, 0, 1); // Must set up to constrained values in bc.

            // Solve for the unknown displacements:::
            Matrix<double> s = Kff.Solve(ff.Subtract(Kfp.Multiply(up)));
            u.SetSubMatrix(nbdof, 0, s); // Set displacements back.

            System.Console.WriteLine(u);
            // Permute back u
            u.PermuteRows(perm);

            System.Console.WriteLine("U after permut");
            System.Console.WriteLine(K);
            System.Console.WriteLine(u);
            System.Console.WriteLine(f);

            // Get reaction forces:
            Matrix<double> Q = K.Multiply(u).Subtract(f);

            ResultFEM result = new ResultFEM(u, Q);

            return result;

        }

        public static double[][][] CoordXtr(int[][] edof, Matrix<double> coord, int[][]dof, int nen)
        {
            int[] coordNr = new int[nen];
            //The structure is [Element nr][x,y,z][node1,node2...]
            double[][][] AllElCoord = new double[edof.Length][][];

            int ndof = dof[0].Length;

            //Loop over all elements
            for (int edofRow = 0; edofRow < edof.Length; edofRow++)
            {

                coordNr[0] = (edof[edofRow][0] - 1) / ndof;
                coordNr[1] = (edof[edofRow][3] - 1) / ndof;                

                //Array with all 
                AllElCoord[edofRow] = new double[coord.ColumnCount][];

                for (int direct = 0; direct < coord.ColumnCount; direct++)
                {
                    AllElCoord[edofRow][direct] = new double[nen];

                    for (int j = 0; j < nen; j++)
                    {
                        AllElCoord[edofRow][direct][j] = coord[coordNr[j], direct];
                    }
                }
            }

            return AllElCoord;
        }

        public static double[][] ExtractElDisp(int[][] edof, double[]a)
        {
            double[][] ed = new double[edof.Length][];

            int i = 0;
            foreach (int[] edofRow in edof)
            {
                ed[i] = new double[edofRow.Length];

                for (int j = 0; j < edofRow.Length; j++)
                {
                    ed[i][j] = a[edofRow[j]-1];
                }
                    
                i++;
            }

            return ed;
        }

        public static double getElementLength(double[][] elCoord)
        {
            double squareSum = 0.0;

            //Loop over x,y,z
            for (int i = 0; i < elCoord.Length; i++)
            {
                squareSum += Math.Pow(elCoord[i][1] - elCoord[i][0],2);
            }

            return Math.Sqrt(squareSum);
        }

        public static void GetGaussPoints(int p, out double[] gaussPoint, out double[] weight)
        {
            double[] newGaussPoint = new double[p + 1];
            double[] newWeight = new double[p + 1];
            if (p == 1)
            {
                newGaussPoint = new double[] { -1.0 / System.Math.Sqrt(3.0), 1.0 / System.Math.Sqrt(3.0) };
                newWeight = new double[] { 1.0, 1.0 };
            }
            else if (p == 2)
            {
                newGaussPoint = new double[] { -System.Math.Sqrt(3.0 / 5.0), 0.0, System.Math.Sqrt(3.0 / 5.0) };
                newWeight = new double[] { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };
            }
            else if (p == 3)
            {
                newGaussPoint = new double[] { -.861136311594053, -.339981043584856, .339981043584856, .861136311594053 };
                newWeight = new double[] { .347854845137454, .652145154862546, .652145154862546, .347854845137454 };
            }
            else if (p == 4)
            {
                newGaussPoint = new double[] { -1.0 / 3.0 * System.Math.Sqrt(5.0 + 2.0 * System.Math.Sqrt(10.0 / 7.0)), -1.0 / 3.0 * System.Math.Sqrt(5.0 - 2.0 * System.Math.Sqrt(10.0 / 7.0)), 0.0, 1.0 / 3.0 * System.Math.Sqrt(5.0 - 2.0 * System.Math.Sqrt(10.0 / 7.0)), 1.0 / 3.0 * System.Math.Sqrt(5.0 + 2.0 * System.Math.Sqrt(10.0 / 7.0)) };
                newWeight = new double[] { (322.0 - 13.0 * System.Math.Sqrt(70.0)) / 900.0, (322.0 + 13.0 * System.Math.Sqrt(70.0)) / 900.0, 128.0 / 225.0, (322.0 + 13.0 * System.Math.Sqrt(70.0)) / 900.0, (322.0 - 13.0 * System.Math.Sqrt(70.0)) / 900.0 };
            }
            else if (p == 5)
            {
                newGaussPoint = new double[] { -.932469514203152, -.661209386466265, -.238619186083197, .238619186083197, .661209386466265, .932469514203152 };
                newWeight = new double[] { .171324492379170, .360761573048139, .467913934572691, .467913934572691, .360761573048139, .171324492379170 };
            }
            gaussPoint = newGaussPoint;
            weight = newWeight;
        }
    
        public static double [] bar3gs(double[][] ec,double[] ep, double[] ed)
        {
            double[] result = new double[2];
            //Computes the Green-Lagrange strain and corresponding normal force in a
            //three dimensional bar element
            //
            // OUTPUT:
            // [es, ee]

            //Initial length squared
            double l02 = Math.Pow(getElementLength(ec),2);

            //Difference in displacement at nodes
            Vector u = new DenseVector(3);

            for (int i = 0; i < 3; i++)
            {
                u[i] = ed[i] - ed[i+3];
            }

            //Bar vector in undeformed configuration
            Vector x0 = new DenseVector(3);

            for (int i = 0; i < 3; i++ )
            {
                //The structure is [x,y,z][node1,node2...]
                x0[i] = ec[i][0] - ec[i][1];
            }

            //Green-lagrange strain
            result[1] = 1 / l02 * (x0.ToRowMatrix().Multiply(u.ToColumnMatrix())[0, 0] + 0.5 * (u.ToRowMatrix().Multiply(u.ToColumnMatrix()))[0, 0]);

            //Normal force
            result[0] = ep[0] * ep[1] * result[1];

            return result;
        }    

        
        public static double[] bar3gfBarnes(double[][] ec,double[] ep, double[] ed, double ten)
        {
            // Compute the internal force vector for a cable element according to:
            // FORM-FINDING  AND  ANALYSIS  OF  PRESTRESSED NETS  AND  MEMBRANES 
            // by M. R. BARNES (1988)
            //ed = [a_1 a_2 ... a_6]; Element nodal displacements
            //ep = [E A]; element properties
            //ten = [ten]; Specified initial element tension

            //OUTPUT:
            //ef = F_int = [f_1 f_2 ... f_6]';
            //ec format: [x,y,z][node1,node2...]
            //ed format: [node1x, node1y, node1z, node2x...]

            double l0 = getElementLength(ec); 
            double[][] ecUpdated = updateEcWithDisp(ec, ed);
            double l = getElementLength(ecUpdated); //Current length

            //Current element force
            double tm = ten + (ep[0] * ep[1] / l0) * (l - l0);

            //Direction vector
            Vector v = new DenseVector(3);

            for (int i = 0; i < 3; i++)
            {
                v[i] = ecUpdated[i][0] - ecUpdated[i][1]; 
            }

            double vl = v.Norm(2);
            Vector x = new DenseVector(6);


            for (int i = 0; i < 3; i++)
            {
                x[i] = v[i] / vl;
            }

            //Reverse direction vector
            for (int i = 0; i < 3; i++)
            {
                x[i+3] = -v[i] / vl;
            }
            

            return x.Multiply(tm).ToArray();

        }

        public static double[][] updateEcWithDisp(double[][] ec, double[] ed)
        {
            //Update ec according do displacements
            double[][] ecUpdated = new double[ec.Length][];
            
            for (int i = 0; i < ec.Length; i++)
            {
                ecUpdated[i] = new double[ec[i].Length];

                for (int j = 0; j < ec[i].Length; j++)
                {
                    ecUpdated[i][j] = ec[i][j] + ed[i + ec.Length * j];
                }
            }
            return ecUpdated;
        }
    }

}
