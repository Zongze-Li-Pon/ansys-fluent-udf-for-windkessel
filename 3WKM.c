// This function needs 10 UDMs. The order of the UDM is:
// 0 -> RSA Qprev, 1 -> RCCA Qprev, 2 -> DAo Qprev, 3 -> LCCA Qprev, 4 -> LSA Qprev,
// 5 -> RSA Pnext, 6 -> RCCA Pnext, 7 -> DAo Pnext, 8 -> LCCA Pnext, 9 -> LSA Pnext.
// This order will be used in this function: F_UDMI(f, t, iWhich) where iWhich gives which UDM is.

// In execute at end macro, ID is necessary for searching the thread, this ID can be found in graphic
// interface of Fluent and may be different for different case, check your ID and input them here.
// The order is RSA, RCCA, DAo, LCCA, LSA.
#include "udf.h"

int threadID[5] = {31, 32, 30, 33, 34};

// Define WKM parameters, order is RSA, RCCA, DAo, LCCA, LSA.
real R1Clinical[5] = {1.907122, 1.854744, 0.262995, 3.497029, 2.839197};
real R2Clinical[5] = {2.362653, 2.448623, 0.347096, 4.661569, 3.72411};
real CClinical[5] = {0.194057, 0.203041, 1.42907, 0.107935, 0.132415};

DEFINE_EXECUTE_AT_END(updateUDMs)
{
    // Get fluid domain.
    Domain* d = Get_Domain(1);
    
    // Variable to be printed out at the end.
    real Qpout[5];
    real Qprevpout[5];
    real Ppout[5];
    real Pnextpout[5];
	// Loop over all openings to find Pnext.
    int i;
    for(i = 0; i < sizeof(threadID)/sizeof(threadID[0]); i++)
    {
        // Get thread with certain ID.
        Thread* t = Lookup_Thread(d, threadID[i]);

        //Define variables used to calculate Pnext.
        face_t f;   // face_t is a type in Fluent, an individual face on boundary.
        real Q, P, area, Qprev, Qg, Pnext, R1, R2, C, k1, k2, k3, k4;    // real is a type in Fluent, double type.
        Q = 0;
        P = 0;
        area = 0;
        Qprev = 0;

        // This macro will loop all faces on this boundary, in this macro, Q and P are got.
        begin_f_loop(f, t)
        {
            real A[3];
            F_AREA(A, f, t);    // This function finds area vector.
            Q += F_FLUX(f, t)/F_R(f, t);    // Mass flowrate/density.
            P += F_P(f, t)*NV_MAG(A);   // NV_MAG() gets magnitude.
            area += NV_MAG(A);  // Get area of this face.
            Qprev += F_UDMI(f, t, i);    // Get previous Q.
        }
        end_f_loop(f, t);
        P /= area;  // Area-weighted pressure.

        // Get Q gradient using first order backward derivative.
        Qg = (Q-Qprev)/CURRENT_TIMESTEP;

        // Find Pnext.
        R1 = R1Clinical[i]*133322368.;    // Convert unit from clinical to SI.
        R2 = R2Clinical[i]*133322368.;    // Convert unit from clinical to SI.
        C = CClinical[i]/133322368.;    // Convert unit from clinical to SI.
        k1 = (R1+R2)*Q/(R2*C)+R1*Qg-P/(R2*C);    // Runge Kutta 4th is used to find dp.
        k2 = (R1+R2)*Q/(R2*C)+R1*Qg-(P+k1*CURRENT_TIMESTEP/2.)/(R2*C);
        k3 = (R1+R2)*Q/(R2*C)+R1*Qg-(P+k2*CURRENT_TIMESTEP/2.)/(R2*C);
        k4 = (R1+R2)*Q/(R2*C)+R1*Qg-(P+k3*CURRENT_TIMESTEP)/(R2*C);
        Pnext = P+1./6.*CURRENT_TIMESTEP*(k1+2.*k2+2.*k3+k4);

        // Update Pnext and Qprev, save them into UDMs.
        begin_f_loop(f, t)
        {
            F_UDMI(f, t, i) = F_FLUX(f, t)/F_R(f, t);
            F_UDMI(f, t, i+5) = Pnext;
        }
        end_f_loop(f, t);

        // Record variables to print out in console.
        Qpout[i] = Q;
        Qprevpout[i] = Qprev;
        Ppout[i] = P;
        Pnextpout[i] = Pnext;
    }

    // Print out current Q and next P.
    printf("Windkessel model updated, current Q, previous Q, P, and Pnext for RSA, RCCA, DAo, LCCA, and LSA are:\n");
    printf("%e, %e, %e, %e, %e\n", Qpout[0], Qpout[1], Qpout[2], Qpout[3], Qpout[4]);
    printf("%e, %e, %e, %e, %e\n", Qprevpout[0], Qprevpout[1], Qprevpout[2], Qprevpout[3], Qprevpout[4]);
    printf("%e, %e, %e, %e, %e\n", Ppout[0], Ppout[1], Ppout[2], Ppout[3], Ppout[4]);
    printf("%e, %e, %e, %e, %e\n", Pnextpout[0], Pnextpout[1], Pnextpout[2], Pnextpout[3], Pnextpout[4]);
}

DEFINE_PROFILE(rsawkm, t, i)
{
	// Apply Pnext onto pressure boundary condition.
    face_t f;
    begin_f_loop(f, t)
    {
        F_PROFILE(f, t, i) = F_UDMI(f, t, 5);
    }
    end_f_loop(f, t);
}
DEFINE_PROFILE(rccawkm, t, i)
{
	// Apply Pnext onto pressure boundary condition.
    face_t f;
    begin_f_loop(f, t)
    {
        F_PROFILE(f, t, i) = F_UDMI(f, t, 6);
    }
    end_f_loop(f, t);
}
DEFINE_PROFILE(daowkm, t, i)
{
	// Apply Pnext onto pressure boundary condition.
    face_t f;
    begin_f_loop(f, t)
    {
        F_PROFILE(f, t, i) = F_UDMI(f, t, 7);
    }
    end_f_loop(f, t);
}
DEFINE_PROFILE(lccawkm, t, i)
{
	// Apply Pnext onto pressure boundary condition.
    face_t f;
    begin_f_loop(f, t)
    {
        F_PROFILE(f, t, i) = F_UDMI(f, t, 8);
    }
    end_f_loop(f, t);
}
DEFINE_PROFILE(lsawkm, t, i)
{
	// Apply Pnext onto pressure boundary condition.
    face_t f;
    begin_f_loop(f, t)
    {
        F_PROFILE(f, t, i) = F_UDMI(f, t, 9);
    }
    end_f_loop(f, t);
}