

#define ASIZE 250
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include "mpi.h"
#include <stdlib.h>
#include <math.h>
#include <err.h>
#include "time.h"

/* X11 data */
static Display *dpy;
static Window win;
static XImage *bitmap;
static char *bitmap_data;
static GC gc;
static Atom wmDeleteMessage;


#define A(x, i, j)  x[(i)+ ((j) * ASIZE)]

static void get_colour(double sv, char *r, char *g, char *b)
{
	int c1 = (sv + 1.0) * 128.0;

	if (c1 > 255) c1 = 255;
	if (c1 < 0) c1 = 0;

	/* Write rgb value */
	*b = (char) 255;
	*g = c1;
	*r = c1;
}

static void pstep(double *restrict a, double *restrict b, double *restrict c, int j)
{
	int i;
	double o, m;

	/* Stable as long as lambda <= 0.5 */
	double lambda = 0.5;

	for (i = 0; i < ASIZE; i++)
	{
		m = A(b, i, j);
		o = A(a, i, j);

		/* Top */
		if (!j)
		{
			if (!i)
			{
				/* Corner */
				A(c, i, j) = 0;
			}
			else if (i == ASIZE - 1)
			{
				/* Corner */
				A(c, i, j) = 0;
			}
			else
			{
				/* Edge */
				A(c, i, j) = lambda * (-m + A(b, i, j + 1) + o - A(a, i, j + 1)
					+ 0.5 * (A(b, i - 1, j) - 2*m + A(b, i + 1, j))) - o + 2*m;
			}
		}

		/* Bottom */
		else if (j == ASIZE - 1)
		{
			if (!i)
			{
				/* Corner */
				A(c, i, j) = 0;
			}
			else if (i == ASIZE - 1)
			{
				/* Corner */
				A(c, i, j) = 0;
			}
			else
			{
				/* Edge */
				A(c, i, j) = lambda * (-m + A(b, i, j - 1) + o - A(a, i, j - 1)
					+ 0.5 * (A(b, i - 1, j) - 2*m + A(b, i + 1, j))) - o + 2*m;
			}
		}

		/* Middle */
		else
		{
			if (!i)
			{
				/* Edge */
				A(c, i, j) = lambda * (-m + A(b, i + 1, j) + o - A(a, i + 1, j)
					+ 0.5 * (A(b, i, j - 1) - 2*m + A(b, i, j + 1))) - o + 2*m;
			}
			else if (i == ASIZE - 1)
			{
				/* Edge */
				A(c, i, j) = lambda * (-m + A(b, i - 1, j) + o - A(a, i - 1, j)
					+ 0.5 * (A(b, i, j - 1) - 2*m + A(b, i, j + 1))) - o + 2*m;
			}
			else
			{
				/* Middle */
				A(c, i, j) = lambda*(A(b, i - 1, j) + A(b, i + 1, j) - 4*m
					+ A(b, i, j - 1) + A(b, i, j + 1)) - o + 2*m;
			}
		}

		get_colour(A(c, i, j), &bitmap->data[i*4 + j*ASIZE*4+2],
				&bitmap->data[i*4 + j*ASIZE*4+1], &bitmap->data[i*4 + j*ASIZE*4]);
	}
}

static void step(double *restrict a, double *restrict b, double *restrict c)
{
	int i, j;
	double o, m;

	/* Stable as long as lambda <= 0.5 */
	double lambda = 0.5;

	for (j = 0; j < ASIZE; j++)
	{
		for (i = 0; i < ASIZE; i++)
		{
			m = A(b, i, j);
			o = A(a, i, j);

			/* Top */
			if (!j)
			{
				if (!i)
				{
					/* Corner */
					A(c, i, j) = 0;
				}
				else if (i == ASIZE - 1)
				{
					/* Corner */
					A(c, i, j) = 0;
				}
				else
				{
					/* Edge */
					A(c, i, j) = lambda * (-m + A(b, i, j + 1) + o - A(a, i, j + 1)
						+ 0.5 * (A(b, i - 1, j) - 2*m + A(b, i + 1, j))) - o + 2*m;
				}
			}

			/* Bottom */
			else if (j == ASIZE - 1)
			{
				if (!i)
				{
					/* Corner */
					A(c, i, j) = 0;
				}
				else if (i == ASIZE - 1)
				{
					/* Corner */
					A(c, i, j) = 0;
				}
				else
				{
					/* Edge */
					A(c, i, j) = lambda * (-m + A(b, i, j - 1) + o - A(a, i, j - 1)
						+ 0.5 * (A(b, i - 1, j) - 2*m + A(b, i + 1, j))) - o + 2*m;
				}
			}

			/* Middle */
			else
			{
				if (!i)
				{
					/* Edge */
					A(c, i, j) = lambda * (-m + A(b, i + 1, j) + o - A(a, i + 1, j)
						+ 0.5 * (A(b, i, j - 1) - 2*m + A(b, i, j + 1))) - o + 2*m;
				}
				else if (i == ASIZE - 1)
				{
					/* Edge */
					A(c, i, j) = lambda * (-m + A(b, i - 1, j) + o - A(a, i - 1, j)
						+ 0.5 * (A(b, i, j - 1) - 2*m + A(b, i, j + 1))) - o + 2*m;
				}
				else
				{
					/* Middle */
					A(c, i, j) = lambda*(A(b, i - 1, j) + A(b, i + 1, j) - 4*m
						+ A(b, i, j - 1) + A(b, i, j + 1)) - o + 2*m;
				}
			}

			get_colour(A(c, i, j), &bitmap->data[i*4 + j*ASIZE*4+2],
				&bitmap->data[i*4 + j*ASIZE*4+1], &bitmap->data[i*4 + j*ASIZE*4]);
		}
	}
}

static void exit_x11(void)
{
	XDestroyImage(bitmap);
	XDestroyWindow(dpy, win);
	XCloseDisplay(dpy);
}

static void init_x11(int size)
{
	unsigned long white, black;

	char name[128] = "Wave Equation";
	char *n = name;
	int depth;

	Visual *visual;

	Status st;

	XTextProperty tp;

	/* Attempt to open the display */
	dpy = XOpenDisplay(NULL);

	/* Failure */
	if (!dpy) exit(0);

   white = WhitePixel(dpy,DefaultScreen(dpy));
   black = BlackPixel(dpy,DefaultScreen(dpy));


	win = XCreateSimpleWindow(dpy, DefaultRootWindow(dpy),
			0, 0, ASIZE, ASIZE, 0, black, white);

	/* We want to be notified when the window appears */
	XSelectInput(dpy, win, StructureNotifyMask);

	/* Make it appear */
	XMapWindow(dpy, win);

	while (1)
	{
		XEvent e;
		XNextEvent(dpy, &e);
		if (e.type == MapNotify) break;
	}

	st = XStringListToTextProperty(&n, 1, &tp);
	if (st) XSetWMName(dpy, win, &tp);

	/* Wait for the MapNotify event */
	XFlush(dpy);

	depth = DefaultDepth(dpy, DefaultScreen(dpy));
	visual = DefaultVisual(dpy, DefaultScreen(dpy));

	if ((depth != 32) && (depth != 24)) errx(1, "Need 32bit colour screen!\n");

	/* Make bitmap */
	bitmap = XCreateImage(dpy, visual, depth,
		ZPixmap, 0, bitmap_data,
		size, size, 32, 0);

	/* Init GC */
	gc = XCreateGC(dpy, win, 0, NULL);
	XSetForeground(dpy, gc, black);

	XSelectInput(dpy, win, ExposureMask | KeyPressMask | StructureNotifyMask);

	wmDeleteMessage = XInternAtom(dpy, "WM_DELETE_WINDOW", False);
	XSetWMProtocols(dpy, win, &wmDeleteMessage, 1);
}

static int test_exit(void)
{
	XEvent event;
	KeySym key;
	char text[255];

	/* Handle pending events */
	while (XPending(dpy))
	{
		/* Get the pending event */
		XNextEvent(dpy, &event);

		/* Press 'q' to quit */
		if ((event.type == KeyPress) &&
			XLookupString(&event.xkey, text, 255, &key, 0) == 1)
		{
			if (text[0] == 'q') return 1;
		}

		/* Or simply close the window */
		if ((event.type == ClientMessage) &&
			((Atom) event.xclient.data.l[0] == wmDeleteMessage))
		{
			return 1;
		}
	}

	return 0;
}

int main(int argc, char **argv)
{
	double *a, *b, *c;

	int rank;
	int nranks;
	int y1, y2;
	int i, j;

	int x, y;

	int count = 0;

	int uprank, downrank;
	MPI_Request up_rq, down_rq;

	int *bitmap_counts;
	int *bitmap_displ;

	int exiting;

	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &nranks);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	y1 = ASIZE * rank / nranks;
	y2 = (ASIZE * (rank + 1) / nranks);

	if (!rank)
	{
		uprank = MPI_PROC_NULL;
		downrank = 1;

		/* Make bitmap */
		bitmap_data = malloc(ASIZE * ASIZE * 4);
		if (!bitmap_data) errx(1, "Couldn't allocate bitmap\n");

		/* Init data for bitmap gather */
		bitmap_counts = malloc(nranks * sizeof(int));
		bitmap_displ = malloc(nranks * sizeof(int));

		bitmap_displ[0] = 0;
		for (i = 0; i < nranks; i++)
		{
			bitmap_counts[i] = 4 * ASIZE * ((ASIZE * (i + 1) / nranks) - ASIZE * i / nranks);
			if (i != nranks - 1) bitmap_displ[i + 1] = bitmap_displ[i] + bitmap_counts[i];
		}

		/* Make a window */
		init_x11(ASIZE);

		/* Init rng */
		srand(time(NULL));
	}
	else
	{
		if (rank != nranks - 1)
		{
			uprank = rank - 1;
			downrank = rank + 1;
		}
		else
		{
			uprank = rank - 1;
			downrank = MPI_PROC_NULL;
		}

		bitmap_data = malloc((y2 - y1) * ASIZE * 4);
		if (!bitmap_data) errx(1, "Couldn't allocate bitmap\n");

		/* Shift bitmap_data */
		bitmap_data -= y1 * ASIZE * 4;

		bitmap_counts = NULL;
		bitmap_displ = NULL;
	}

	/* Allocate shifted arrays using void * arithmetic */
	a = malloc(sizeof(double) * (y2 - y1 + 2) * ASIZE);
	b = malloc(sizeof(double) * (y2 - y1 + 2) * ASIZE);
	c = malloc(sizeof(double) * (y2 - y1 + 2) * ASIZE);

	if (!a || !b || !b) errx(1, "Failed to allocate arrays\n");

	/* Shift the arrays so that we can use zero-based offsets */
	a -= (y1 - 1) * ASIZE;
	b -= (y1 - 1) * ASIZE;
	c -= (y1 - 1) * ASIZE;

	/* Init my part of the array + overlap */
	for (j = y1 - 1; j <= y2; j++)
	{
		for (i = 0; i < ASIZE; i++)
		{
			double dx = i - ASIZE/2;
			double dy = j - ASIZE/2;

			double d2 = (dx*dx+dy*dy) / 20;

			A(a, i, j) = 10 * exp(-d2);
			A(b, i, j) = A(a, i, j);
		}
	}

	while (1)
	{
		/* Calculate boundaries as early as possible, and then send them */
		pstep(a, b, c, y1);
		MPI_Isend(&A(c, 0, y1), ASIZE, MPI_DOUBLE, uprank, 0, MPI_COMM_WORLD, &up_rq);

		pstep(a, b, c, y2 - 1);
		MPI_Isend(&A(c, 0, y2 - 1), ASIZE, MPI_DOUBLE, downrank, 0, MPI_COMM_WORLD, &down_rq);

		for (j = y1 + 1; j < y2 - 1; j++)
		{
			pstep(a, b, c, j);
		}

		/* Wait for the results */
		MPI_Recv(&A(c, 0, y1 - 1), ASIZE, MPI_DOUBLE, uprank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&A(c, 0, y2), ASIZE, MPI_DOUBLE, downrank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Wait(&up_rq, MPI_STATUS_IGNORE);
		MPI_Wait(&down_rq, MPI_STATUS_IGNORE);

		/* Calculate boundaries as early as possible, and then send them */
		pstep(b, c, a, y1);
		MPI_Isend(&A(a, 0, y1), ASIZE, MPI_DOUBLE, uprank, 0, MPI_COMM_WORLD, &up_rq);

		pstep(b, c, a, y2 - 1);
		MPI_Isend(&A(a, 0, y2 - 1), ASIZE, MPI_DOUBLE, downrank, 0, MPI_COMM_WORLD, &down_rq);

		for (j = y1 + 1; j < y2 - 1; j++)
		{
			pstep(b, c, a, j);
		}

		/* Wait for the results */
		MPI_Recv(&A(a, 0, y1 - 1), ASIZE, MPI_DOUBLE, uprank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&A(a, 0, y2), ASIZE, MPI_DOUBLE, downrank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Wait(&up_rq, MPI_STATUS_IGNORE);
		MPI_Wait(&down_rq, MPI_STATUS_IGNORE);

		/* Calculate boundaries as early as possible, and then send them */
		pstep(c, a, b, y1);
		MPI_Isend(&A(b, 0, y1), ASIZE, MPI_DOUBLE, uprank, 0, MPI_COMM_WORLD, &up_rq);

		pstep(c, a, b, y2 - 1);
		MPI_Isend(&A(b, 0, y2 - 1), ASIZE, MPI_DOUBLE, downrank, 0, MPI_COMM_WORLD, &down_rq);

		for (j = y1 + 1; j < y2 - 1; j++)
		{
			pstep(c, a, b, j);
		}

		/* Wait for the results */
		MPI_Recv(&A(b, 0, y1 - 1), ASIZE, MPI_DOUBLE, uprank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&A(b, 0, y2), ASIZE, MPI_DOUBLE, downrank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		MPI_Wait(&up_rq, MPI_STATUS_IGNORE);
		MPI_Wait(&down_rq, MPI_STATUS_IGNORE);

		/* Send image information to root */
		MPI_Gatherv(&bitmap_data[y1*ASIZE*4], ASIZE*4*(y2 - y1), MPI_BYTE,
					bitmap_data, bitmap_counts, bitmap_displ, MPI_BYTE, 0, MPI_COMM_WORLD);

		if (!rank)
		{
			/* Display it */
			XPutImage(dpy, win, gc, bitmap,
						0, 0, 0, 0, ASIZE, ASIZE);

			XFlush(dpy);

			exiting = test_exit();
		}
		MPI_Bcast(&exiting, 1, MPI_INT, 0, MPI_COMM_WORLD);
		if (exiting) break;

		count++;

		/* Add another drop? */
		if (count == 100)
		{
			count = 0;

			/* Decide where */
			if (!rank)
			{
				x = (((ASIZE - 101.0) * rand()) / RAND_MAX) + 50.0;
				y = (((ASIZE - 101.0) * rand()) / RAND_MAX) + 50.0;
			}

			/* Tell other ranks */
			MPI_Bcast(&x, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast(&y, 1, MPI_INT, 0, MPI_COMM_WORLD);

			for (j = y1 - 1; j <= y2; j++)
			{
				for (i = 0; i < ASIZE; i++)
				{
					double dx = i - x;
					double dy = j - y;

					double d2 = (dx*dx+dy*dy) / 20;

					A(a, i, j) += 10 * exp(-d2);
					A(b, i, j) += 10 * exp(-d2);
				}
			}
		}
	}

	/* Invert shifts */
	a += (y1 - 1) * ASIZE;
	b += (y1 - 1) * ASIZE;
	c += (y1 - 1) * ASIZE;
	bitmap_data += y1 * ASIZE * 4;

	free(a);
	free(b);
	free(c);

	/* Done! */
	if (!rank)
	{
		/* X will free our bitmap for us */
		exit_x11();
	}
	else
	{
		free(bitmap_data);
	}

	MPI_Finalize();

	return 0;
}
