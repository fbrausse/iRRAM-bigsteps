
#define _POSIX_C_SOURCE	200809L	/* getline() */
#define _XOPEN_SOURCE	700	/* M_PI */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>		/* size_t, etc. */
#include <math.h>		/* sin()/cos() */
#include <cairo.h>
#include <SDL.h>
#include <sys/time.h>		/* gettimeofday() */

#define PROG_TITLE	"Pendulum Visualizer"
#define INPUT_BUFFER_SZ	(1 << 20)
#define PENDULUM_LENGTH	200	/* px */
#define PENDULUM_RADIUS	30	/* px */

struct line_buf {
	char *line;
	size_t sz;
};

static int read_record(
	FILE *f, struct line_buf *buf,
	unsigned n_fields, const unsigned *field_ids,
	double *values
) {
	ssize_t len;

loop:
	while ((len = getline(&buf->line, &buf->sz, f)) > 0) {
		const char *tok;
		unsigned id;
		unsigned valid_fields = 0;
		for (id = 1, tok = strtok(buf->line, " \t");
		     tok && valid_fields < n_fields;
		     id++, tok = strtok(NULL, " \t")) {
			if (tok[0] == '#')	/* remaining line is comment */
				goto loop;
			for (unsigned id_idx = 0; id_idx < n_fields; id_idx++) {
				if (field_ids[id_idx] != id)
					continue;
				if (sscanf(tok, "%lf", values + id_idx) != 1) {
					fprintf(stderr,
						"cannot interpret %u:'%s' as "
						"double; ignoring line\n",
						id, tok);
					goto loop;
				}
				valid_fields++;
			}
		}
		if (valid_fields == n_fields)
			return 0;
		fprintf(stderr, "not enough fields on line, ignoring\n");
	}

	return 1;
}

struct pt {
	double x, y;
};

static void paint_bobbel(cairo_t *cr, const struct pt *d, double r)
{
	cairo_pattern_t *pat;
	pat = cairo_pattern_create_radial(
			d->x - .3 * r, d->y - .3 * r, 0.1 * r,
			d->x - .3 * r, d->y - .3 * r, 1.3 * r);
	cairo_pattern_add_color_stop_rgba(pat, 0, 1, 1, 1, 1);
	cairo_pattern_add_color_stop_rgba(pat, 1, 0, 0, 0, 1);

	cairo_save(cr);
	cairo_set_source(cr, pat);
	cairo_arc(cr, d->x, d->y, r, 0, 2 * M_PI);
	cairo_fill(cr);
	cairo_restore(cr);

	cairo_pattern_destroy(pat);
}

static void paint_pendulums(
	cairo_t *cr, const struct pt *c, const double *v, unsigned n,
	double r, double len
) {
	if (n == 0)
		return;

	struct pt d = {
		c->x + sin(*v) * len,
		c->y + cos(*v) * len,
	};

	/* paint line c -> d */
	cairo_save(cr);
	cairo_move_to(cr, c->x, c->y);
	cairo_line_to(cr, d.x, d.y);
	cairo_stroke(cr);
	cairo_restore(cr);

	paint_pendulums(cr, &d, v+1, n-1, r, len);

	paint_bobbel(cr, &d, r);
}

struct frame {
	unsigned bpp, w, h, s;
};

int main(int argc, char **argv)
{
	struct line_buf buf;
	unsigned n_fields = 0;
	unsigned i;
	const char *tok;
	char *field_spec;
	double time_scale = 1.0;

	/* adjust input buffer to attenuate delays due to iRRAM reiterations */
	static char input_buffer[INPUT_BUFFER_SZ];
	setvbuf(stdin, input_buffer, _IOFBF, sizeof(input_buffer));

	/* option parsing */
	double display_rate = INFINITY;
	int opt;
	while ((opt = getopt(argc, argv, ":hr:s:")) != -1)
		switch (opt) {
		case 'r': display_rate = atof(optarg); break;
		case 's': time_scale   = atof(optarg); break;
		case 'h': goto usage;
		case ':':
			fprintf(stderr, "option '-%c' needs an argument\n",
				optopt);
			exit(1);
		case '?':
			fprintf(stderr, "unknown option '-%c'\n", optopt);
			exit(1);
		}
	if (optind != argc-1) {
usage:
		fprintf(stderr,
			"usage: %s [-s <time-scale>] [-r <display-rate>] "
			"<t:phi1[:phi2[:...]]>\n",
			argv[0]);
		exit(1);
	}
	field_spec = argv[optind];

	/* parse field_spec */
	for (tok = field_spec-1; tok; tok = strchr(tok + 1, ':'))
		n_fields++;
	unsigned field_ids[n_fields];
	for (i=0, tok = strtok(field_spec, ":"); tok;
	     i++, tok = strtok(NULL, ":"))
		field_ids[i] = atoi(tok);

	/* init video output */
	struct frame f;
	f.bpp = 32;
	f.h = f.w = 2 * ((n_fields-1) * PENDULUM_LENGTH + PENDULUM_RADIUS) + 20;

	if (SDL_Init(SDL_INIT_VIDEO) < 0) {
		fprintf(stderr, "error to init'ing SDL: %s\n", SDL_GetError());
		exit(1);
	}
	SDL_Surface *sdl_sf = SDL_SetVideoMode(f.w, f.h, f.bpp,
	                                       SDL_SWSURFACE | SDL_DOUBLEBUF);
	if (!sdl_sf) {
		fprintf(stderr, "error init'ing SDL video surface: %s\n",
			SDL_GetError());
		exit(1);
	}
	f.s = sdl_sf->pitch;

	SDL_WM_SetCaption(PROG_TITLE " - press ESC to quit", PROG_TITLE);

	double input_last_t = 0;

	struct timeval tv;
	gettimeofday(&tv, NULL);
	double exec_last_t = tv.tv_sec + tv.tv_usec / 1e6;

	double values[n_fields];
	memset(&buf, 0, sizeof(buf));
	for (unsigned rec=0;
	     !read_record(stdin, &buf, n_fields, field_ids, values);
	     rec++) {
		double t = values[0] * time_scale;
		double input_delta_t = t - input_last_t;
		if (input_delta_t < 1.0/display_rate)
			continue;

		gettimeofday(&tv, NULL);
		double exec_t = tv.tv_sec + tv.tv_usec / 1e6;
		double exec_delta_t = exec_t - exec_last_t;

		if (input_delta_t > exec_delta_t)
			usleep((input_delta_t - exec_delta_t) * 1e6);

		SDL_Event event;
		while (SDL_PollEvent(&event))
			switch (event.type) {
			case SDL_QUIT:
				goto done;
			case SDL_KEYDOWN:
				switch (event.key.keysym.sym) {
				case SDLK_ESCAPE:
				case SDLK_q:
					goto done;
				default:
					break;
				}
				break;
			default:
				break;
			}

		SDL_LockSurface(sdl_sf);
		cairo_surface_t *sf = cairo_image_surface_create_for_data(
				sdl_sf->pixels,
				CAIRO_FORMAT_ARGB32, f.w, f.h, f.s);
		cairo_t *cr = cairo_create(sf);
		cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
		cairo_paint(cr);
		cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);

		paint_pendulums(cr, &(struct pt){ f.w/2.0, f.h/2.0 },
		                values+1, n_fields-1,
		                PENDULUM_RADIUS, PENDULUM_LENGTH);

		/* draw some record infos on the frame */
		cairo_move_to(cr, 10, f.h - 5);
		static char legend_buf[64];
		snprintf(legend_buf, sizeof(legend_buf)-1,
			"Step %u, time: %.6e", rec, values[0]);
		legend_buf[sizeof(legend_buf)-1] = '\0';
		cairo_show_text(cr, legend_buf);
		for (unsigned field = 1; field < n_fields; field++) {
			snprintf(legend_buf, sizeof(legend_buf)-1,
				", phi%u: %.e", field, values[field]);
			legend_buf[sizeof(legend_buf)-1] = '\0';
			cairo_show_text(cr, legend_buf);
		}

		cairo_destroy(cr);
		cairo_surface_destroy(sf);
		SDL_UnlockSurface(sdl_sf);
		SDL_Flip(sdl_sf);		/* update display */

		input_last_t = t;
		gettimeofday(&tv, NULL);
		exec_last_t = tv.tv_sec + tv.tv_usec / 1e6;
	}
done:
	free(buf.line);

	SDL_Quit();

	return 0;
}
