#include "embed.h"
#include "navier-stokes/centered.h"
#include "navier-stokes/perfs.h"

scalar f[];

#ifdef PECLET
#include "tracer.h"
scalar *tracers = {f};
#include "diffusion.h"
#else // !PECLET
#include "vof.h"
scalar *interfaces = {f};
#endif

#include "view.h"

const double Re = 0.00164;     // horizontal Reynolds
const double Uv = 1;       // relative vertical velocity
const double muv = 10. ; // relative viscosity of the vertical fluid
const double Hv = 1.;       // relative vertical channel diameter
#ifdef PECLET
const double Pe = PECLET;
#endif

#define poiseuille(x) (6. * ((x) + 0.5) * (0.5 - (x)))

u.n[left] = dirichlet(cs[] * poiseuille(y));
p[left] = neumann(0.);
pf[left] = neumann(0.);
f[left] = (cs[] > 0.);

u.n[right] = neumann(0.);
p[right] = dirichlet(0.);
pf[right] = dirichlet(0.);

u.n[bottom] = dirichlet(+Uv * cs[] * poiseuille(x));
p[bottom] = neumann(0.);
pf[bottom] = neumann(0.);
f[bottom] = 0.;

u.n[top] = dirichlet(-Uv * cs[] * poiseuille(x));
p[top] = neumann(0.);
pf[top] = neumann(0.);
f[top] = 0.;

u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);

const int maxlevel = 9;

face vector mut[];

int main()
{
    L0 = 3.000001; // just to make sure that the solid boundaries do not exactly coincide with cell boundaries
    DT = HUGE;
    origin(-1.25, -L0 / 2.);
    N = 128;
    mu = mut;
#ifdef PECLET
    f.gradient = minmod2;
#endif

#if 0 // c'est du bricolage à améliorer 
  TOLERANCE = HUGE;
  NITERMIN = 10;
  NITERMAX = 10;
#endif

    run();
}

#ifdef PECLET
face vector CD[];
#endif

event properties(i++)
{
    double muh = 1. / Re, muV = muv * muh;
    foreach_face()
    {
        mut.x[] = fm.x[] * ((muh - muV) * ((clamp(f[], 0, 1) + clamp(f[-1], 0, 1)) / 2.) + muV);
#ifdef PECLET
        CD.x[] = fm.x[] / Pe;
#endif
    }
}

event init(t = 0)
{
    solid(cs, fs, union(intersection(-y + 1. / 2., +y + 1. / 2.), intersection(-x + Hv / 2., +x + Hv / 2.)));
    foreach ()
        f[] = (cs[] > 0.) * (x < 0.5) * (fabs(y) - Delta / 2. < 0.5);
}

#ifdef PECLET
event tracer_diffusion(i++)
{
    diffusion(f, dt, CD);
}
#endif // PECLET

event adapt(t = 2. / Uv; i++)
{
    const double mmax = 1e-2, fmax = 1e-2, uemax = 0.1;
    adapt_wavelet({cs, f, u}, (double[]){mmax, fmax, uemax, uemax}, maxlevel);
}

event thread_width(i++)
{
    double s = 0., sf = 0.;
    foreach (reduction(+ : s) reduction(+ : sf))
        if (x > L0 / 2.)
            s += dv(), sf += f[] * dv();
    fprintf(stderr, "%g %g\n", t, sf / s);
}

event movie(t = 4. / Uv)
{
    view(tx = -0.118, ty = -0.000, tz = -3.038,
         width = 918, height = 890);

    draw_vof(c = "f", lc = {1, 1, 1}, lw = 2);
    draw_vof(c = "cs", s = "fs", edges = true, filled = -1, fc = {1, 1, 1});
    squares(color = "sqrt(u.x*u.x + u.y*u.y)", spread = -1, linear = true);
    box();
    //  vectors (u = "u", scale = 0.001);
    save("u.png");

    draw_vof(c = "f", lc = {1, 1, 1}, lw = 2);
    draw_vof(c = "cs", s = "fs", edges = true, filled = -1, fc = {1, 1, 1});
    squares(color = "sqrt(u.x*u.x + u.y*u.y)", spread = -1, linear = true);
    box();
    cells();
    save("cells.png");

    output_facets(f, stdout);

    FILE *fp = fopen("profile", "w");
    for (double y = -0.5; y < 0.51; y += 0.01)
        fprintf(fp, "%g %g\n", y, interpolate(u.x, L0 / 2., y));
    fclose(fp);
}