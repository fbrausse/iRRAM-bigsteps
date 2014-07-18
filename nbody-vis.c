
/* x86_64-pc-mingw32-gcc -static -s -std=c11 -O2 -Wall -DNDEBUG -Wextra -pedantic -Wno-unused-parameter -I/usr/x86_64-pc-mingw32/usr/include/{SDL,cairo} -D_GNU_SOURCE=1 -o nbody.exe nbody-vis.c ring-buf.c -mwindows -lmingw32 -lSDLmain -lSDL -ldxguid -lcairo -lws2_32 -lfreetype -lpixman-1 -lz -lwinmm */

#define _XOPEN_SOURCE		500		/* M_PI, usleep() */

#include <stdlib.h>
#include <string.h>
#include <unistd.h>				/* size_t, etc. */
#include <math.h>				/* sin()/cos() */
#include <cairo.h>
#include <SDL.h>
#include <sys/time.h>				/* gettimeofday() */
#include <fcntl.h>				/* fcntl() */
#include <errno.h>				/* errno, EAGAIN */
#include <assert.h>
#if defined(_WIN64) || defined(_WIN32)
# include <winsock2.h>				/* select() et al. */
#endif

#include "ring-buf.h"

#define TRAIL_QUANTIZED		1
#define TRAIL_MESHED		2

#define PROG_TITLE		"n-Body Visualizer"
#define INPUT_BUFFER_RESERVE	(1 << 22)	/* byte */
#define BODY_RADIUS		(1.0/6.0)	/* input distance unit */
#define TRAIL_WIDTH		(BODY_RADIUS/4)	/* input distance unit */
#define CENTER_X		0		/* input distance unit */
#define CENTER_Y		0		/* input distance unit */
#define WIDTH			800		/* px */
#define HEIGHT			800		/* px */
#define SCALE			72		/* px/(input distance unit) */
#define VIEW_SCROLL_RATIO	(1.0/6.0)
#define VIEW_SCALE_RATIO	(2.0/3.0)
#define TIME_SCALE_RATIO	(2.0/3.0)
#define EVENT_LOOP_SLEEP	500		/* micro seconds */
//#define TRAIL_QUANT_STEPS	256
//#define TRAIL_METHOD		TRAIL_MESHED
#define TRAIL_MAX_SEGMENTS	2000
#define TRAIL_MAX_DIST		4.5		/* px */
#define GAMMA			2.2
#define AUTO_ZOOM_FRAME_BORDER	20
#define LIGHT_COLOR_R		255
#define LIGHT_COLOR_G		255
#define LIGHT_COLOR_B		255
#define COLOR_MAX		255

#define FRAME_CAIRO_FORMAT	CAIRO_FORMAT_RGB24
#define FRAME_CAIRO_ASHIFT	24
#define FRAME_CAIRO_RSHIFT	16
#define FRAME_CAIRO_GSHIFT	8
#define FRAME_CAIRO_BSHIFT	0
#define FRAME_CAIRO_AMASK	(0x00 << FRAME_CAIRO_ASHIFT)
#define FRAME_CAIRO_RMASK	(0xff << FRAME_CAIRO_RSHIFT)
#define FRAME_CAIRO_GMASK	(0xff << FRAME_CAIRO_GSHIFT)
#define FRAME_CAIRO_BMASK	(0xff << FRAME_CAIRO_BSHIFT)

/* check the above settings and define some macros used by the code later on */

#ifndef TRAIL_METHOD
# define TRAIL_METHOD		TRAIL_QUANTIZED
#endif

#if TRAIL_METHOD == TRAIL_QUANTIZED
# ifndef TRAIL_QUANT
#  ifdef TRAIL_QUANT_STEPS
#   define TRAIL_QUANT(v)	(int)((v) * (TRAIL_QUANT_STEPS - 1.0) + .5)
#  else
#   define TRAIL_QUANT(v)	(v)
#  endif
# endif
#elif TRAIL_METHOD == TRAIL_MESHED
# if !(CAIRO_VERSION_MAJOR == 1 && CAIRO_VERSION_MINOR >= 12)
#  error TRAIL_METHOD == TRAIL_MESHED requires cairo >= 1.12
# endif
#else
# error unknown TRAIL_METHOD, supported are TRAIL_QUANTIZED, TRAIL_MESHED
#endif

#ifndef FRAME_GRAY
# if LIGHT_COLOR_R == LIGHT_COLOR_G && LIGHT_COLOR_G == LIGHT_COLOR_B
#  define FRAME_GRAY		1
# else
#  define FRAME_GRAY		0
# endif
#endif

#ifndef M_PI
# define M_PI			atan2(+0.0, -0.0)
#endif

static int parse_record(
	const char *line, unsigned n_fields, const unsigned *field_ids,
	double *values
) {
	unsigned id;
	unsigned valid_fields = 0;
	int a, n;

	for (id = 1; *line && valid_fields < n_fields; id++) {
		unsigned id_idx;
		line += strspn(line, " \t");
		if (line[0] == '#')	/* remaining line is comment */
			return 1;
		a = strcspn(line, " \t");
		for (id_idx = 0; id_idx < n_fields; id_idx++) {
			if (field_ids[id_idx] != id)
				continue;
			if (sscanf(line, "%lf%n", values + id_idx, &n) < 1 ||
			    n != a) {
				fprintf(stderr,
					"cannot interpret %u:'%.*s' as "
					"double; ignoring line\n",
					id, (int)a, line);
				return 1;
			}
			valid_fields++;
		}
		line += a;
	}
	if (valid_fields >= n_fields)
		return 0;
	fprintf(stderr, "not enough fields on line, ignoring\n");
	return 1;
}

static int read_record(
	int fd, struct ring_buf *buf, struct ring_buf *line,
	unsigned n_fields, const unsigned *field_ids,
	double *values
) {
	size_t newline = 0;

	if (rb_read(buf, fd) == -1 && errno != EAGAIN)
		return -1;

	while (1) {
		ssize_t rd;
		char *q = rb_chr(buf, '\n');
		if (q) {
#if 1
			newline = rb_ptrdiff(buf, q);
#else
			size_t p = q - buf->buf;
			newline = p >= buf->idx ? p - buf->idx
			                        : buf->sz + p - buf->idx;
#endif
		}

		if (newline) {
#if 1
			rb_reset(line);
#else
			line->idx = line->valid = 0;
#endif
			rb_ensure_sz(line, newline+1);
			rb_tfer(line, buf, newline);
			rb_putc(line, '\0');
			rb_getc(buf); /* eat '\n' (if available) */
			/* fprintf(stderr, "'%s'\n", line->buf); */
			if (!parse_record(line->buf, n_fields, field_ids, values))
				return 1;
			newline = 0;
			if (buf->valid)
				continue;
		}

		if (buf->valid == buf->sz)
			rb_ensure_sz(buf, buf->sz < INPUT_BUFFER_RESERVE
			                  ? INPUT_BUFFER_RESERVE
			                  : 2 * buf->sz);

		rd = rb_read(buf, fd);/*
		fprintf(stderr,
			"buf after read %zd: idx: %zu, valid: %zu, sz: %zu, "
			"newline: %zu\n",
			rd, buf->idx, buf->valid, buf->sz, newline);*/
		if (rd < 0)
			return errno == EAGAIN ? 0 : -1;
		if (rd == 0) {
			if (!buf->valid)
				return -1;
			newline = buf->valid;
		} else {
			newline = 0;
		}
	}
}

struct pt {
	double x, y;
};

static inline struct pt * pt_add(struct pt *a, const struct pt *b)
{
	a->x += b->x;
	a->y += b->y;
	return a;
}

static inline struct pt * pt_neg(struct pt *a)
{
	a->x = -a->x;
	a->y = -a->y;
	return a;
}

static inline struct pt * pt_sub(struct pt *a, const struct pt *b)
{
	a->x -= b->x;
	a->y -= b->y;
	return a;
}

static inline struct pt * pt_muld(struct pt *a, double s)
{
	a->x *= s;
	a->y *= s;
	return a;
}

static inline double pt_norm_square(const struct pt *a)
{
	return a->x * a->x + a->y * a->y;
}

static inline double pt_norm(const struct pt *a)
{
	return hypot(a->x, a->y);
}

static inline struct pt * pt_min(struct pt *min, const struct pt *b)
{
	if (b->x < min->x) min->x = b->x;
	if (b->y < min->y) min->y = b->y;
	return min;
}

static inline struct pt * pt_max(struct pt *max, const struct pt *b)
{
	if (b->x > max->x) max->x = b->x;
	if (b->y > max->y) max->y = b->y;
	return max;
}

/* returns determinant d of [1 ax ay; 1 bx by; 1 cx cy]:
 * d < 0 => abc cw, d > 0 => abc ccw, d == 0 => abc collinear */
static double pt_orientation(
	const struct pt *a, const struct pt *b, const struct pt *c
) {
	return (b->x - a->x)*(c->y - a->y) - (c->x - a->x)*(b->y - a->y);
}

static void paint_bobbel(cairo_t *cr, const struct pt *d, double r)
{
	cairo_pattern_t *pat;
	pat = cairo_pattern_create_radial(
			d->x - .3 * r, d->y - .3 * r, 0.1 * r,
			d->x - .3 * r, d->y - .3 * r, 1.3 * r);
	cairo_pattern_add_color_stop_rgb(pat, 0,
	                                 (double)LIGHT_COLOR_R / COLOR_MAX,
	                                 (double)LIGHT_COLOR_G / COLOR_MAX,
	                                 (double)LIGHT_COLOR_B / COLOR_MAX);
	cairo_pattern_add_color_stop_rgb(pat, 1, 0, 0, 0);

	cairo_save(cr);
	cairo_set_source(cr, pat);
	cairo_arc(cr, d->x, d->y, r, 0, 2 * M_PI);
	cairo_fill(cr);
	cairo_restore(cr);

	cairo_pattern_destroy(pat);
}

__attribute__((format(printf,2,3)))
static int cr_printf(cairo_t *cr, const char *fmt, ...)
{
	va_list ap;
	int n;
	va_start(ap, fmt);
	n = vsnprintf(NULL, 0, fmt, ap);
	va_end(ap);
	if (n < 0)
		return n;

	char buf[n+1];
	va_start(ap, fmt);
	vsnprintf(buf, n+1, fmt, ap);
	va_end(ap);

	cairo_show_text(cr, buf);
	return n;
}

struct frame {
	unsigned w, h;
	cairo_surface_t *sf;
	SDL_Surface *sdl_sf;
	int sdl_flags;
};

struct pos_it {
	struct pos_it *next;
	double t;
	struct pt pts[];
};

struct pos_queue {
	struct pos_it *head;
	struct pos_it *tail;
	size_t n;
};

static void pos_queue_update(
	struct pos_queue *q, unsigned n_fields, const double *values,
	double trail_fade_t
) {
	struct pos_it *it = NULL;
	unsigned i;
	while (q->head && values[0] - q->head->t > trail_fade_t) {
		free(it);
		it = q->head;
		q->head = it->next;
		if (!q->head)
			q->tail = NULL;
		q->n--;
	}
	if (!it)
		it = malloc(sizeof(struct pos_it) +
		            sizeof(struct pt) * (n_fields/2));
	it->next = NULL;
	it->t = values[0];
	for (i=0; i<n_fields/2; i++) {
		it->pts[i].x = values[1+i*2+0];
		it->pts[i].y = values[1+i*2+1];
	}
	if (q->tail)
		q->tail->next = it;
	else
		q->head = it;
	q->tail = it;
	q->n++;
}

/* cw */
static struct pt orth(const struct pt *s, const struct pt *t, double scale)
{
	struct pt r = { -(t->y - s->y), (t->x - s->x) };
	scale /= hypot(r.x, r.y);
	r.x *= scale;
	r.y *= scale;
	return r;
}

#define CR_PATH_OBJS_CLOSE	(1+0)
#define CR_PATH_OBJS_MOVE_TO	(1+1)
#define CR_PATH_OBJS_LINE_TO	(1+1)
#define CR_PATH_OBJS_CURVE_TO	(1+3)

static void cr_path_close(cairo_path_data_t *d)
{
	d->header.type = CAIRO_PATH_CLOSE_PATH;
	d->header.length = 1;
}

static void cr_path_move_to(cairo_path_data_t *d, const struct pt *p)
{
	d->header.type = CAIRO_PATH_MOVE_TO;
	d->header.length = 2;
	d++;
	d->point.x = p->x;
	d->point.y = p->y;
}

static void cr_path_line_to(cairo_path_data_t *d, const struct pt *p)
{
	d->header.type = CAIRO_PATH_LINE_TO;
	d->header.length = 2;
	d++;
	d->point.x = p->x;
	d->point.y = p->y;
}

#if TRAIL_METHOD == TRAIL_QUANTIZED
/* FIXME: endpoints of adjoining line segments get more ink due to painting
 *        those twice */
static int paint_trail(
	cairo_t *cr,
	double t, unsigned n_pts, const struct pos_queue *q, double trail_fade_t
) {
	cairo_save(cr);
	cairo_push_group(cr);
	cairo_set_operator(cr, CAIRO_OPERATOR_SOURCE);
	cairo_set_line_width(cr, TRAIL_WIDTH);
	/*cairo_set_line_cap(cr, CAIRO_LINE_CAP_ROUND);*/
	cairo_set_line_join(cr, CAIRO_LINE_JOIN_BEVEL);
	const struct pos_it *prev, *it = q->head;
	unsigned iterations = 0;
	while ((prev = it) && (it = it->next)) {
		iterations++;
		double trail_fraction = (t - it->t) / trail_fade_t;
		double c = 1.0 - pow(trail_fraction, 1.0/GAMMA);
		unsigned n = 1, i;
		double d = c;
		while (it->next) {
			double next_trail_fraction = (t - it->next->t) / trail_fade_t;
			d = 1.0 - pow(next_trail_fraction, 1.0/GAMMA);
			if (TRAIL_QUANT(c) != TRAIL_QUANT(d))
				break;
			it = it->next;
			n++;
		}
#if 0
		for (i=0; i<n_pts; i++) {
			cairo_pattern_t *pat = cairo_pattern_create_linear(
					prev->pts[i].x,
					prev->pts[i].y,
					(it->next ? it->next : it)->pts[i].x,
					(it->next ? it->next : it)->pts[i].y);
			cairo_pattern_add_color_stop_rgb(pat, 0.0, c, c, c);
			cairo_pattern_add_color_stop_rgb(pat, 1.0, d, d, d);
			cairo_set_source(cr, pat);
			cairo_pattern_destroy(pat);
			cairo_move_to(cr, prev->pts[i].x, prev->pts[i].y);
			const struct pos_it *it2 = prev->next;
			for (unsigned j=0; j<=n && it2; j++, it2 = it2->next)
				cairo_line_to(cr, it2->pts[i].x, it2->pts[i].y);
			cairo_stroke(cr);
		}
#elif 1
		cairo_set_source_rgba(cr,0,0,0,c);
		for (i=0; i<n_pts; i++) {
			const struct pos_it *it2 = prev->next;
			unsigned j;
			cairo_move_to(cr, prev->pts[i].x, prev->pts[i].y);
			for (j=0; j<=n && it2; j++, it2 = it2->next)
				cairo_line_to(cr, it2->pts[i].x, it2->pts[i].y);
		}
		cairo_stroke(cr);
#else
		cairo_set_source_rgba(cr,0,0,0,c);
		for (i=0; i<n_pts; i++) {
			cairo_path_data_t path[n*4*(1+1)];
			unsigned i0 = 0, i1 = n*4*(1+1);
			unsigned j;

			const struct pos_it *it0 = prev;
			const struct pos_it *it1 = it0->next;
			const struct pos_it *it2 = it1->next;
			struct pt o = orth(it0->pts+i, it1->pts+i, TRAIL_WIDTH/2.0);
			struct pt p = orth(it1->pts+i, it2->pts+i, TRAIL_WIDTH/2.0);
			struct pt q = o;
			pt_muld(pt_add(&q, &p), .5);
			struct pt p1 = it0->pts[i];
			struct pt p2 = it1->pts[i];
			struct pt p3 = it1->pts[i];
			struct pt p4 = it0->pts[i];
			pt_add(&p1, &o);
			pt_sub(&p4, &o)

			// cairo_move_to(cr, p1.x, p1.y);
			cr_path_move_to(&path[i0], &p1), i0 += CR_PATH_OBJS_MOVE_TO;
			cr_path_close(&path[i1 -= CR_PATH_OBJS_CLOSE]);
			cr_path_line_to(&path[i1 -= CR_PATH_OBJS_LINE_TO], &p4);

			double orient = pt_orientation(it0->pts+i, it1->pts+i, it2->pts+i);
			if (orient > 0.0) { /* ccw */
			} else if (orient < -0.0) { /* cw */
			} else { /* collinear */
			}
		}
		cairo_fill(cr);
#endif
	}
	cairo_pop_group_to_source(cr);
	cairo_paint(cr);
	cairo_restore(cr);
	return iterations;
}
#elif TRAIL_METHOD == TRAIL_MESHED
static void cr_mesh_patch(
	cairo_pattern_t *mesh,
	const struct pt *side0_center, const struct pt *side0_off, double a0,
	const struct pt *side2_center, const struct pt *side2_off, double a2
) {
	const struct pt *s = side0_center;
	const struct pt *t = side2_center;
	const struct pt *o = side0_off;
	const struct pt *p = side2_off;

	cairo_mesh_pattern_begin_patch(mesh);
	cairo_mesh_pattern_move_to(mesh, s->x + o->x, s->y + o->y);
	cairo_mesh_pattern_line_to(mesh, s->x - o->x, s->y - o->y);
	cairo_mesh_pattern_line_to(mesh, t->x - p->x, t->y - p->y);
	cairo_mesh_pattern_line_to(mesh, t->x + p->x, t->y + p->y);
	cairo_mesh_pattern_set_corner_color_rgba(mesh, 0, 0, 0, 0, a0);
	cairo_mesh_pattern_set_corner_color_rgba(mesh, 1, 0, 0, 0, a0);
	cairo_mesh_pattern_set_corner_color_rgba(mesh, 2, 0, 0, 0, a2);
	cairo_mesh_pattern_set_corner_color_rgba(mesh, 3, 0, 0, 0, a2);
	cairo_mesh_pattern_end_patch(mesh);
}

/* FIXME: intersections between parts of the same path are painted
 *        non-antialiased, because mesh patches have non-AA borders */
static int paint_trail(
	cairo_t *cr,
	double t, unsigned n_pts, const struct pos_queue *q, double trail_fade_t
) {
	const struct pos_it *prev, *it, *next;
	unsigned iterations = 0, i;
	cairo_save(cr);
	cairo_set_line_width(cr, TRAIL_WIDTH);
	/*cairo_set_line_cap(cr, CAIRO_LINE_CAP_ROUND);*/
	cairo_set_line_join(cr, CAIRO_LINE_JOIN_BEVEL);
	for (i=0; i<n_pts; i++) {
		it = q->head;
		if (!it || !(next = it->next))
			break;
		cairo_pattern_t *mesh = cairo_pattern_create_mesh();
		cairo_move_to(cr, it->pts[i].x, it->pts[i].y);

		while ((prev = it, it = it->next) != NULL) {
			double col0 = 1.0 - pow((t - prev->t) / trail_fade_t, 1.0/GAMMA);
			double col1 = 1.0 - pow((t - it->t) / trail_fade_t, 1.0/GAMMA);
			const struct pt *s = prev->pts + i;
			const struct pt *t = it->pts + i;

			cairo_line_to(cr, t->x, t->y);

			/* orth on (s -> t) with length TRAIL_WIDTH/2.0 */
			struct pt o = orth(s, t, TRAIL_WIDTH / 2.0);
			/* define a col0-to-col1 gradient patch that completely
			 * covers the current line segment */
			cr_mesh_patch(mesh, s, &o, col0, t, &o, col1);
			if (!it->next)
				continue;
			/* orth on (t -> u) with length TRAIL_WIDTH/2.0 */
			struct pt p = orth(t, it->next->pts+i, TRAIL_WIDTH/2.0);
			/* define a solid col1 patch that completely covers the
			 * (bevel) join of the current line with the next one */
			cr_mesh_patch(mesh, t, &o, col1, t, &p, col1);
		}

		cairo_set_source(cr, mesh);
		cairo_stroke(cr);
		cairo_pattern_destroy(mesh);
	}
	cairo_restore(cr);
	return iterations;
}
#endif

#include <inttypes.h>

static void paint_frame(
	const struct frame *f, const cairo_matrix_t *pt2view,
	unsigned n_fields, const double *values, unsigned rec,
	const struct pos_queue *q, double trail_fade_t, double rb_filled,
	double exec_delta_t, int nrec
) {
	unsigned n_pts = n_fields/2, i;
	struct timeval tv, tw;
	double x, y, w, h;
	int trail_n = 0;

	gettimeofday(&tv, NULL);

	cairo_t *cr = cairo_create(f->sf);
	cairo_set_operator(cr, CAIRO_OPERATOR_SOURCE);
	cairo_set_source_rgba(cr, 1.0, 1.0, 1.0, 1.0);
	cairo_paint(cr);
	cairo_set_operator(cr, CAIRO_OPERATOR_OVER);

	cairo_set_matrix(cr, pt2view);

	if (q)
		trail_n = paint_trail(cr, values[0], n_pts, q, trail_fade_t);

	for (i=0; i<n_pts; i++) {
		x = values[1+2*i+0];
		y = values[1+2*i+1];
		paint_bobbel(cr, &(struct pt){ x, y }, BODY_RADIUS);
	}
	x = f->w / 2.0;
	y = f->h / 2.0;
	w = f->w / 2.0;
	h = f->h / 2.0;
	cairo_device_to_user(cr, &x, &y);
	cairo_device_to_user_distance(cr, &w, &h);

	cairo_identity_matrix(cr);

	/* draw some record infos on the frame */
	cairo_set_source_rgb(cr, 0, 0, 0);
	cairo_move_to(cr, 5, f->h - 5);
	cr_printf(cr, "Step %u, time: %.6e", rec, values[0]);
	for (i=0; i<n_pts; i++)
		cr_printf(cr, ", x%u: %.6e y%u: %.6e",
			i, values[1+2*i+0], i, values[1+2*i+1]);

	cairo_move_to(cr, 5, 13);
	cr_printf(cr, "view: [%.6e +/- %.6e] x [%.6e +/- %.6e]", x,w, y,h);

	gettimeofday(&tw, NULL);
	cr_printf(cr,
		", frame creation: %02.1f ms, loop: %02.1f ms, queue: %" PRIuMAX ":%u"
		", buf: %4.1f %%, nrec: %d",
		(tw.tv_sec - tv.tv_sec) * 1e3 + (tw.tv_usec - tv.tv_usec) / 1e3,
		exec_delta_t * 1e3,
		(uintmax_t)q->n, trail_n, rb_filled * 100, nrec);

	cairo_destroy(cr);
	cairo_surface_flush(f->sf);
}

static void frame_set_surface(struct frame *f)
{
	cairo_format_t format = FRAME_CAIRO_FORMAT;
	if (f->sf) {
		if ((unsigned)cairo_image_surface_get_width(f->sf) == f->w &&
		    (unsigned)cairo_image_surface_get_height(f->sf) == f->h &&
		    cairo_image_surface_get_format(f->sf) == format)
			return;
		cairo_surface_destroy(f->sf);
	}
	f->sf = cairo_image_surface_create(format, f->w, f->h);
}

static int sdl_set_video_mode(struct frame *f)
{
	if (!SDL_SetVideoMode(f->w, f->h, 0, f->sdl_flags))
		return 1;
	if (f->sdl_sf)
		SDL_FreeSurface(f->sdl_sf);
	frame_set_surface(f);
	/* these SDL surfaces never need locking */
	f->sdl_sf = SDL_CreateRGBSurfaceFrom(
			cairo_image_surface_get_data(f->sf),
			f->w, f->h, 32, cairo_image_surface_get_stride(f->sf),
			FRAME_CAIRO_RMASK, FRAME_CAIRO_GMASK, FRAME_CAIRO_BMASK,
			FRAME_CAIRO_AMASK);
	return 0;
}

static void sdl_update_display(const struct frame *f)
{
	SDL_Surface *screen = SDL_GetVideoSurface();
	SDL_BlitSurface(f->sdl_sf, NULL, screen, NULL);
	SDL_Flip(screen);
}

struct at {
	double x0;
	double s;
};

static double affine_transform(const struct at *a, double p)
{
	return a->x0 + p * a->s;
}

static double affine_inverse_transform(const struct at *a, double q)
{
	return (q - a->x0) / a->s;
}

#define EV_QUIT		(1 << 0)
#define EV_REPAINT	(1 << 1)

static int process_event(
	const SDL_Event *ev, struct frame *f, cairo_matrix_t *a, double t,
	struct at *time_scale, int *auto_zoom
) {
	struct pt tmp = { 0, 0 };
	cairo_matrix_t inv_a;
	double v0;
	memcpy(&inv_a, a, sizeof(inv_a));
	cairo_matrix_invert(&inv_a);

	switch (ev->type) {
	case SDL_QUIT:
		return EV_QUIT;
	case SDL_KEYDOWN:
		switch (ev->key.keysym.sym) {
		case SDLK_ESCAPE:
		case SDLK_q:
			return EV_QUIT;
		case SDLK_a:
			*auto_zoom = !*auto_zoom;
			return EV_REPAINT;
		case SDLK_PLUS:
		case SDLK_KP_PLUS:
			if (!isinf(t)) {
				v0 = affine_inverse_transform(time_scale, t);
				time_scale->s *= TIME_SCALE_RATIO;
				time_scale->x0 = t - time_scale->s * v0;
			} else {
				time_scale->s *= TIME_SCALE_RATIO;
			}
			break;
		case SDLK_MINUS:
		case SDLK_KP_MINUS:
			if (!isinf(t)) {
				v0 = affine_inverse_transform(time_scale, t);
				time_scale->s /= TIME_SCALE_RATIO;
				time_scale->x0 = t - time_scale->s * v0;
			} else {
				time_scale->s /= TIME_SCALE_RATIO;
			}
			break;
		case SDLK_0:
		case SDLK_KP0:
			if (!isinf(t)) {
				v0 = affine_inverse_transform(time_scale, t);
				time_scale->s = 1;
				time_scale->x0 = t - time_scale->s * v0;
			} else {
				time_scale->s = 1;
			}
			break;
		case SDLK_LEFT:
			tmp.x = f->w * VIEW_SCROLL_RATIO;
			cairo_matrix_transform_distance(&inv_a, &tmp.x, &tmp.y);
			cairo_matrix_translate(a, tmp.x, tmp.y);
			*auto_zoom = 0;
			return EV_REPAINT;
		case SDLK_RIGHT:
			tmp.x = f->w * -VIEW_SCROLL_RATIO;
			cairo_matrix_transform_distance(&inv_a, &tmp.x, &tmp.y);
			cairo_matrix_translate(a, tmp.x, tmp.y);
			*auto_zoom = 0;
			return EV_REPAINT;
		case SDLK_UP:
			tmp.y = f->h * VIEW_SCROLL_RATIO;
			cairo_matrix_transform_distance(&inv_a, &tmp.x, &tmp.y);
			cairo_matrix_translate(a, tmp.x, tmp.y);
			*auto_zoom = 0;
			return EV_REPAINT;
		case SDLK_DOWN:
			tmp.y = f->h * -VIEW_SCROLL_RATIO;
			cairo_matrix_transform_distance(&inv_a, &tmp.x, &tmp.y);
			cairo_matrix_translate(a, tmp.x, tmp.y);
			*auto_zoom = 0;
			return EV_REPAINT;
		default:
			break;
		}
		break;
	case SDL_MOUSEMOTION:
		if (ev->motion.state & SDL_BUTTON(1)) {
			tmp.x = ev->motion.xrel;
			tmp.y = ev->motion.yrel;
			cairo_matrix_transform_distance(&inv_a, &tmp.x, &tmp.y);
			cairo_matrix_translate(a, tmp.x, tmp.y);
			*auto_zoom = 0;
			return EV_REPAINT;
		}
		break;
	case SDL_MOUSEBUTTONDOWN:
		switch (ev->button.button) {
		case 4: /* wheel up */
			tmp.x = ev->button.x;
			tmp.y = ev->button.y;
			cairo_matrix_transform_point(&inv_a, &tmp.x, &tmp.y);
			cairo_matrix_translate(a, tmp.x, tmp.y);
			cairo_matrix_scale(a, 1.0/VIEW_SCALE_RATIO,
					   1.0/VIEW_SCALE_RATIO);
			cairo_matrix_translate(a, -tmp.x, -tmp.y);
			*auto_zoom = 0;
			return EV_REPAINT;
		case 5: /* wheel down */
			tmp.x = ev->button.x;
			tmp.y = ev->button.y;
			cairo_matrix_transform_point(&inv_a,
						     &tmp.x, &tmp.y);
			cairo_matrix_translate(a, tmp.x, tmp.y);
			cairo_matrix_scale(a, VIEW_SCALE_RATIO,
					   VIEW_SCALE_RATIO);
			cairo_matrix_translate(a, -tmp.x, -tmp.y);
			*auto_zoom = 0;
			return EV_REPAINT;
		}
		break;
	case SDL_VIDEORESIZE:
		f->w = ev->resize.w;
		f->h = ev->resize.h;
		if (sdl_set_video_mode(f)) {
			fprintf(stderr,
				"error resizing surface to %dx%d: %s\n",
				f->w, f->h, SDL_GetError());
			exit(1);
		}
		return EV_REPAINT;
	default:
		break;
	}

	return 0;
}

static struct timeval * tv_norm(struct timeval *a)
{
	ldiv_t n = ldiv(a->tv_usec, 1000000);
	if (n.rem < 0) {
		n.quot--;
		n.rem += 1000000;
	}
	a->tv_sec += n.quot;
	a->tv_usec = n.rem;
	return a;
}

#if defined(_WIN64) || defined(_WIN32)
# define suseconds_t	signed long
#endif

static inline struct timeval * tv_addi(struct timeval *a, time_t sec, suseconds_t usec)
{
	a->tv_sec += sec;
	a->tv_usec += usec;
	return tv_norm(a);
}

static inline struct timeval * tv_add(struct timeval *a, const struct timeval *b)
{
	return tv_addi(a, b->tv_sec, b->tv_usec);
}

static inline struct timeval * tv_subi(struct timeval *a, time_t sec, suseconds_t usec)
{
	a->tv_sec -= sec;
	a->tv_usec -= usec;
	return tv_norm(a);
}

static inline struct timeval * tv_sub(struct timeval *a, const struct timeval *b)
{
	return tv_subi(a, b->tv_sec, b->tv_usec);
}

static inline int tv_test(const struct timeval *a)
{
	if (a->tv_sec < 0) return -1;
	if (a->tv_sec > 0) return +1;
	if (a->tv_usec < 0) return -1;
	if (a->tv_usec > 0) return +1;
	return 0;
}

static int tv_rem(struct timeval *limit)
{
	struct timeval now;
	gettimeofday(&now, NULL);
	return tv_test(tv_sub(limit, &now));
}

/* calc n-1 finite difference values; store n values in fd
 * 0: point
 * 1: speed
 * 2: acceleration
 * ... */
static void fdiff_update(double *fd, unsigned n, double y, double inv_dt)
{
	unsigned i;

	for (i=0; i<n; i++) {
		double x = fd[i];
		fd[i] = i ? (fd[i-1] - y) * inv_dt : y;
		y = x;
	}
	/*
	struct fdiff2 r;
	r.x   = y;
	r.dx  = (r.x  - fd2->x ) * inv_dt;
	r.ddx = (r.dx - fd2->dx) * inv_dt;
	*fd2  = r;*/
	/*
	fd->s = *pt_muld(pt_sub(&d, &fd->b), inv_dt);
	fd->b = *c;*/
	/*
	fd->s.x = (c->x - fd->b.x) / dt;
	fd->s.y = (c->y - fd->b.y) / dt;
	fd->a = fd->b;
	fd->b = *c;
	fd->dt = dt;
	fd->ds_square = fd->s.x * fd->s.x + fd->s.y * fd->s.y;*/
}

struct zoom_control_state {
	double x0[3], y0[3], x1[3], y1[3];	/* fdiffs: required bbox */
	double scale[3];			/* derived from bbox, needed? */
};

struct bbox {
	struct pt min, max;
};

#define BBOX_INIT	{ { NAN, NAN }, { NAN, NAN } }

/* returns whether bbox a entirely contains bbox b */
static int bbox_contains(const struct bbox *a, const struct bbox *b)
{
	if (b->min.x < a->min.x) return 0;
	if (b->min.y < a->min.y) return 0;
	if (b->max.x > a->max.x) return 0;
	if (b->max.y > a->max.y) return 0;
	return 1;
}

static inline struct pt * bbox_dim(struct pt *ret, const struct bbox *b)
{
	*ret = b->max;
	return pt_sub(ret, &b->min);
}
#if 0
static double bbox_scale(
	const struct pt *a_dim, const struct pt *b_dim
) {
	double sx = a_dim->x / b_dim->x;
	double sy = a_dim->y / b_dim->y;
	return sx < sy ? sx : sy;
}
#endif
static inline void bbox_span_pt(struct bbox *b, const struct pt *p, int empty)
{
	if (empty) {
		memcpy(&b->min, p, sizeof(b->min));
		memcpy(&b->max, p, sizeof(b->max));
	} else {
		pt_min(&b->min, p);
		pt_max(&b->max, p);
	}
}

static inline void bbox_span(struct bbox *a, const struct bbox *b, int empty)
{
	if (empty) {
		memcpy(a, b, sizeof(*a));
	} else {
		pt_min(&a->min, &b->min);
		pt_max(&a->max, &b->max);
	}
}

static inline int bbox_is_empty(const struct bbox *b)
{
	return !(b->min.x <= b->max.x && b->min.y <= b->max.y);
}

static void scene_bbox(
	struct bbox *bbox_device, const cairo_matrix_t *a, const double *values,
	unsigned n, double body_r, const struct pos_it *it, double trail_r
) {
	unsigned i, j;

	if (!n)
		return;

	/* find bbox for center points of all objects in device space */
	for (i=0; i<n; i++) {
		struct pt p = { values[i+i+0], values[i+i+1] };
		cairo_matrix_transform_point(a, &p.x, &p.y);
		bbox_span_pt(bbox_device, &p, !i);
	}

	/* adjust for the radius; this is in user space, so a transformation of
	 * the device unit vectors to user space is done to determine the
	 * appropriate scaling factors for r regarding horizontal and vertical
	 * boundaries;
	 * are these only valid for the current scaling factor of the matrix, so
	 * body_r_device needs to be rescaled accordingly later??? */
	struct pt hori = { 1, 0 }, vert = { 0, 1 };
	cairo_matrix_t inv_a = *a;
	cairo_matrix_invert(&inv_a);
	cairo_matrix_transform_distance(&inv_a, &hori.x, &hori.y);
	cairo_matrix_transform_distance(&inv_a, &vert.x, &vert.y);
	double inv_hori_norm = 1.0 / pt_norm(&hori);
	double inv_vert_norm = 1.0 / pt_norm(&vert);

	struct pt body_r_device = {
		.x = body_r * inv_hori_norm,
		.y = body_r * inv_vert_norm,
	};
	pt_sub(&bbox_device->min, &body_r_device);
	pt_add(&bbox_device->max, &body_r_device);

	if (!it)
		return;

	struct bbox trail_bbox_device;
	for (j=0; it; it = it->next)
		for (i=0; i<n; i++, j++) {
			struct pt p = it->pts[i];
			cairo_matrix_transform_point(a, &p.x, &p.y);
			bbox_span_pt(&trail_bbox_device, &p, !j);
		}
	struct pt trail_r_device = {
		.x = trail_r * inv_hori_norm,
		.y = trail_r * inv_vert_norm,
	};
	pt_sub(&trail_bbox_device.min, &trail_r_device);
	pt_add(&trail_bbox_device.max, &trail_r_device);
	bbox_span(bbox_device, &trail_bbox_device, !j);
}

static void center_view2(
	const cairo_matrix_t *a, const struct bbox *bbox_device,
	struct zoom_control_state *zctl, unsigned w, unsigned h,
	cairo_matrix_t *b
) {
	struct pt bbox_device_dim;
	bbox_dim(&bbox_device_dim, bbox_device);

#if 0
	struct bbox bbox_last = {
		{ zctl->x0[0], zctl->y0[0] },
		{ zctl->x1[0], zctl->y1[0] },
	};
	struct pt bbox_last_dim;
	bbox_dim(&bbox_last_dim, &bbox_last);

	if (bbox_contains(&bbox_last, bbox_device) &&
	    bbox_scale(&bbox_last_dim, &bbox_device_dim) < 1.5
	) {
		bbox_device_min = (struct pt){ zctl->x0[0], zctl->y0[0] };
		bbox_device_max = (struct pt){ zctl->x1[0], zctl->y1[0] };
		bbox_device_dim = bbox_device_max;
	}

	bbox_dim(&bbox_device_dim, bbox_device);
#endif
	struct pt trans = { 0, 0 };
	pt_sub(&trans, &bbox_device->min);

	/* center view */
	if (bbox_device_dim.x < w)
		trans.x += (w - bbox_device_dim.x) / 2.0;
	if (bbox_device_dim.y < h)
		trans.y += (h - bbox_device_dim.y) / 2.0;

	double sx = w / bbox_device_dim.x;
	double sy = h / bbox_device_dim.y;
	double s = sx < sy ? sx : sy;

	cairo_matrix_init_identity(b);
	cairo_matrix_scale(b, s, s);
	cairo_matrix_translate(b, trans.x, trans.y);
}

static void center_view3(
	cairo_matrix_t *pt2view, const double *values, unsigned n,
	struct zoom_control_state *zctl, double dt, const struct pos_it *it,
	const struct frame *f
) {
	struct bbox bbox_device;
	cairo_matrix_t b;
	int fb = AUTO_ZOOM_FRAME_BORDER;

	scene_bbox(&bbox_device, pt2view, values, n, BODY_RADIUS,
	           it, TRAIL_WIDTH / 2.0);

	if (!n)
		return;

	/* update the bounding box for the current time step of length dt */
	double inv_dt = 1.0 / dt;
	fdiff_update(zctl->x0, 3, bbox_device.min.x, inv_dt);
	fdiff_update(zctl->y0, 3, bbox_device.min.y, inv_dt);
	fdiff_update(zctl->x1, 3, bbox_device.max.x, inv_dt);
	fdiff_update(zctl->y1, 3, bbox_device.max.y, inv_dt);

	center_view2(pt2view, &bbox_device, zctl, f->w-2*fb, f->h-2*fb, &b);

	cairo_matrix_translate(&b, fb, fb);
	cairo_matrix_multiply(pt2view, pt2view, &b);
}

static void process_interactively(
	int fd, unsigned n_fields, const unsigned *field_ids, struct frame *f,
	cairo_matrix_t *pt2view, struct at *time_scale,
	double trail_fade_t, double display_rate, int auto_zoom
) {
	struct zoom_control_state zctl;
	struct ring_buf buf;
	struct ring_buf rb;

	/* init video output */
	if (SDL_Init(SDL_INIT_VIDEO | SDL_INIT_NOPARACHUTE) < 0) {
		fprintf(stderr, "error to init'ing SDL: %s\n", SDL_GetError());
		exit(1);
	}
	if (sdl_set_video_mode(f)) {
		fprintf(stderr, "error init'ing SDL video surface: %s\n",
			SDL_GetError());
		exit(1);
	}

	SDL_WM_SetCaption(PROG_TITLE " - press ESC to quit", PROG_TITLE);

	double input_last_t = 0;

	struct timeval tv;
	gettimeofday(&tv, NULL);
	double exec_last_t = tv.tv_sec + tv.tv_usec / 1e6;

	double values[n_fields];
	memset(&buf, 0, sizeof(buf));
	memset(&rb, 0, sizeof(rb));

	double t = -INFINITY;
	int need_new_record = 1;
	unsigned rec = 0;
	struct pos_queue q;
	q.head = NULL;
	q.tail = NULL;
	q.n = 0;
	int delay = 0;
	int nrec = 0;
	while (1) {
		SDL_Event ev;
		while (SDL_PollEvent(&ev))
			if (process_event(&ev, f, pt2view, t, time_scale, &auto_zoom) & EV_QUIT)
				goto done;
		switch (read_record(fd, &rb, &buf, n_fields, field_ids, values)) {
		case -1: goto done;
		case  0: continue;
		case  1: goto event_loop;
		}
	}
event_loop:
	while (1) {
		SDL_Event ev;
		int ev_mask = 0;
		while (SDL_PollEvent(&ev))
			ev_mask |= process_event(&ev, f, pt2view, t,
						 time_scale, &auto_zoom);
		gettimeofday(&tv, NULL);
		double exec_t  = tv.tv_sec + tv.tv_usec / 1e6;
		double exec_delta_t    = exec_t - exec_last_t;
		if (ev_mask & EV_QUIT)
			goto done;
		if (ev_mask & EV_REPAINT) {
			if (auto_zoom)
				center_view3(pt2view, values+1, n_fields/2,
				             &zctl, 0, q.head, f); /* TODO: dt */
			paint_frame(f, pt2view, n_fields, values, rec, &q,
			            trail_fade_t, (double)rb.valid / rb.sz,
			            exec_delta_t, nrec);
			sdl_update_display(f);
			continue;
		}

#if 1
		t = affine_transform(time_scale, values[0]);
		double input_delta_t   = t - input_last_t;
		double display_delta_t = 1.0 / display_rate;

		/* nrec: have next record available
		 * nrec               delay, req new record, update display
		 * 
		 *  no, di < de < dd:   1           1             0        
		 *  no, de < di < dd:   1           1             0        
		 *  no, di < dd < de:   0           1             1!       
		 *  no, de < dd < di:   1           0             0        
		 *  no, dd < di < de:   0           1             1        
		 *  no, dd < de < di:   1           0             0        
		 * yes, di < de < dd:   1           1             0        
		 * yes, de < di < dd:   1           1             0        
		 * yes, di < dd < de:   0           1             0!       
		 * yes, de < dd < di:   1           0             0        
		 * yes, dd < di < de:   0           1             1        
		 * yes, dd < de < di:   1           0             0        
		 *
		 * update display => new record, no delay
		 * de < di        => no update display, no new record
		 * de < dd        => no update display, delay
		 */
#undef min
#undef max
#define min(a,b)	((a) < (b) ? (a) : (b))
#define max(a,b)	((a) > (b) ? (a) : (b))

		if (input_delta_t <= exec_delta_t &&
		    display_delta_t <= (nrec ? input_delta_t : exec_delta_t)) {
			/* update display */
			input_last_t = t;
			exec_last_t = exec_t;
			/* redraw image */
			if (!q.tail ||
			    TRAIL_MAX_SEGMENTS * (values[0] - q.tail->t) >= trail_fade_t)
				pos_queue_update(&q, n_fields, values, trail_fade_t);
			if (auto_zoom)
				center_view3(pt2view, values+1, n_fields/2,
				             &zctl, 0, q.head, f); /* TODO: dt: exec_delta_t? */
			paint_frame(f, pt2view, n_fields, values, rec, &q,
			            trail_fade_t, (double)rb.valid / rb.sz,
			            exec_delta_t, nrec);
			sdl_update_display(f);
			need_new_record = 1;
			delay = 0;
		} else {
			need_new_record = input_delta_t < max(display_delta_t,
			                                      exec_delta_t);
			delay           = exec_delta_t < max(display_delta_t,
			                                     input_delta_t);
		}

		struct timeval delay_until = tv;
		tv_addi(&delay_until, 0, delay ? EVENT_LOOP_SLEEP : 0);
		if (need_new_record) {
			double tmp_values[n_fields];
			nrec = 1;
			while (1) {
				int ret;
				ret = read_record(fd, &rb, &buf, n_fields,
				                  field_ids, tmp_values);
				if (ret == -1)
					goto done;
				if (ret == 1) {
					rec++;
					memcpy(values, tmp_values, sizeof(values));
					break;
				}
				fd_set fds;
				FD_ZERO(&fds);
				FD_SET(fd, &fds);
				struct timeval timeout = delay_until;
				gettimeofday(&tv, NULL);
				tv_sub(&timeout, &tv);
				if (tv_test(&timeout) <= 0 || 
				    !select(fd+1, &fds, NULL, NULL, &timeout)) {
					nrec = 0;
					break;
				}
			}
		}
		if (delay) {
			/*rb_read(&rb, fd);*/
#if 0
			if (tv_rem(&delay_until) > 0)
				usleep(delay_until.tv_sec * 1000000 + delay_until.tv_usec);
#elif 0
			gettimeofday(&tv, NULL);
			tv_sub(&delay_until, &tv);
			// tv_subi(&delay_until, 0, 100);
			if (tv_test(&delay_until) > 0)
				usleep(delay_until.tv_sec * 1000000 + delay_until.tv_usec);
#else
			double sleep_t = display_delta_t - exec_delta_t;
			if (sleep_t * 1e6 > 2 * EVENT_LOOP_SLEEP)
				usleep(EVENT_LOOP_SLEEP);
#endif
		}
		continue;

#else
		if (need_new_record) {
			double tmp_values[n_fields];
			while (1) {
				int ret;
				ret = read_record(fd, &rb, &buf, n_fields,
				                  field_ids, tmp_values);
				if (ret == -1)
					goto done;
				if (ret == 1)
					break;
				fd_set fds;
				FD_ZERO(&fds);
				FD_SET(fd, &fds);
				struct timeval timeout = {0, EVENT_LOOP_SLEEP};
				if (!select(fd+1, &fds, NULL, NULL, &timeout))
					goto event_loop;
			}
			rec++;
			memcpy(values, tmp_values, sizeof(values));
		}

		need_new_record = 1;
		t = affine_transform(&time_scale, values[0]);
		double input_delta_t = t - input_last_t;

		gettimeofday(&tv, NULL);
		double exec_t = tv.tv_sec + tv.tv_usec / 1e6;
		double exec_delta_t = exec_t - exec_last_t;
#if 1
		if (input_delta_t < 1.0/display_rate)
			continue;
#else
		double max_delta_t = exec_delta_t > input_delta_t
		                   ? exec_delta_t : input_delta_t;
		if (max_delta_t < 1.0/display_rate)
			continue;/*
		if (exec_delta_t > input_delta_t) {
			exec_t = t;
			exec_delta_t = input_delta_t;
		}*/
#endif
		if (input_delta_t <= exec_delta_t) {
			input_last_t = t;
			exec_last_t = exec_t;
			/* redraw image */
			pos_queue_update(&q, n_fields, values, trail_fade_t);
			paint_frame(&f, &pt2view, n_fields, values, rec, &q,
			            trail_fade_t, (double)rb.valid / rb.sz);
			emit_frame(&f, t);
			exec_last_t = gettime();
			continue;
		}
		need_new_record = 0;

		double sleep_t = input_delta_t - exec_delta_t;
		if (sleep_t * 1e6 > 2 * EVENT_LOOP_SLEEP)
			usleep(EVENT_LOOP_SLEEP);
#endif
	}
done:
	free(buf.buf);
	free(rb.buf);

	SDL_Quit();
}

enum output_mode {
	SDL,
	PPM,
	Y4M,
};

union output_data {
	enum output_mode mode;
	struct sdl_data {
		enum output_mode mode;
		SDL_Surface *sdl_sf;
		int sdl_flags;
	} sdl;
	struct ppm_data {
		enum output_mode mode;
		FILE *fofd;
		unsigned char *img;
	} ppm;
	struct y4m_data {
		enum output_mode mode;
		FILE *fofd;
		unsigned char *img;
	} y4m;
};
/*
struct output_class {
	void (*begin_stream)(FILE *fofd, union output_data *odat, const frame *f);
	void (*write_frame)(FILE *fofd, union output_data *odat, cairo_surface_t *sf);
	void (*end_stream)(FILE *fofd, union output_data *odat);
};*/

static void rgb2yuv(
	unsigned char *y, unsigned char *u, unsigned char *v,
	unsigned r, unsigned g, unsigned b
) {
#if 1
#define CSCALE(n,bit,fraction) \
	(((uintmax_t)(n) << (bit)) * (fraction) / 255.0 + ((fraction) < 0 ? -.5 : +.5))
/* 8.24 bit fixed point */
#define CSCALE3(min,max,f1,v1,f2,v2,f3,v3) \
	(( \
		(int32_t)CSCALE((max)-(min),24,(f1)) * (v1) + \
		(int32_t)CSCALE((max)-(min),24,(f2)) * (v2) + \
		(int32_t)CSCALE((max)-(min),24,(f3)) * (v3) + \
		((int32_t)1 << 23) \
	) >> 24)

	/* ITU-R-B.601 coefficients */
	*y =  16 + CSCALE3(16,235,+.299   ,r,+.587   ,g,+.114   ,b);
	*u = 128 + CSCALE3(16,240,-.168736,r,-.331264,g,+.5     ,b);
	*v = 128 + CSCALE3(16,240,+.5     ,r,-.418688,g,-.081312,b);
#else
	/* the coefficients of r,g,b below are
	 * 256**3*219/255*{ITU-R-B.601-coeff} */
	*y = 16U + ((uint32_t)(
			4308192U * r +
			8457888U * g +
			1642588U * b +
			(1U << 23)
		) >> 24);
	*u = 128 + ((int32_t)(
			-2431261 * r +
			-4773073 * g +
			+7204334 * b +
			(1 << 23)
		) >> 24);
	*v = 128 + ((int32_t)(
			+7204334 * r +
			-6032736 * g +
			-1171598 * b +
			(1 << 23)
		) >> 24);
#endif
}

static void process_noninteractively(
	int ifd, FILE *fofd, const enum output_mode omode, unsigned n_fields,
	const unsigned *field_ids, struct frame *f, cairo_matrix_t *pt2view,
	struct at *time_scale, double trail_fade_t, double display_rate,
	int auto_zoom
) {
	struct ring_buf rb = RING_BUF_INIT;
	struct ring_buf scratch = RING_BUF_INIT;
	double values[2][n_fields];

	unsigned long rec, frame_nr = 0;
	unsigned s, x, y;
	size_t img_sz;
	const unsigned char *data;
	unsigned char *img, *p, *u, *v;
	fd_set fds;
	struct pos_queue q;
	struct zoom_control_state zctl;
	q.head = NULL;
	q.tail = NULL;
	q.n = 0;

	frame_set_surface(f);

	switch (omode) {
	case PPM:
		img_sz = f->h * f->w * 3;
		if (!(img = malloc(img_sz))) {
			perror("malloc");
			exit(1);
		}
		break;
	case Y4M:
		img_sz = f->h * f->w * 3 / 2;
		if (!(img = malloc(img_sz))) {
			perror("malloc");
			exit(1);
		}
		fprintf(fofd, "YUV4MPEG2 W%u H%u C420jpeg Ip F%u:1 A1:1\n",
			f->w, f->h, (unsigned)round(display_rate));
#if FRAME_GRAY == 1
		memset(img + f->h * f->w, 0x80, f->h * f->w / 2);
#endif
		break;
	}

	while (1) {
		switch (read_record(ifd, &rb, &scratch, n_fields, field_ids,
		                    values[0])) {
		case -1: goto done;
		case  0:
			FD_ZERO(&fds);
			FD_SET(ifd, &fds);
			if (select(ifd+1, &fds, NULL, NULL, NULL) == -1) {
				perror("select");
				goto done;
			}
			continue;
		case  1: break;
		}
		break;
	}
	for (rec=1;;) {
		double *cur_values = values[~rec&1];
		double *next_values = values[rec&1];
		switch (read_record(ifd, &rb, &scratch, n_fields, field_ids,
		                    next_values)) {
		case -1: goto done;
		case  0:
			FD_ZERO(&fds);
			FD_SET(ifd, &fds);
			if (select(ifd+1, &fds, NULL, NULL, NULL) == -1) {
				perror("select");
				goto done;
			}
			continue;
		case  1:
			rec++;
			break;
		}

		double frame_t = frame_nr / display_rate;

		double next_t = affine_transform(time_scale, next_values[0]);
		if (next_t < frame_t) {
			if (!q.tail) {
				pos_queue_update(&q, n_fields, next_values, trail_fade_t);
			} else {
				unsigned i;
				for (i=0; i<n_fields/2; i++) {
					struct pt p = q.tail->pts[i];
					struct pt r = {
						next_values[1+i+i+0],
						next_values[1+i+i+1],
					};
					pt_sub(&r, &p);
					cairo_matrix_transform_distance(pt2view, &r.x, &r.y);
					if (pt_norm_square(&r) >= TRAIL_MAX_DIST * TRAIL_MAX_DIST) {
						pos_queue_update(&q, n_fields, next_values, trail_fade_t);
						break;
					}
				}
			}
			continue;
		}

		double cur_t = affine_transform(time_scale, cur_values[0]);
		const double *nearest_values;
		if (fabs(frame_t - cur_t) < fabs(next_t - frame_t)) {
			nearest_values = cur_values;
		} else {
			nearest_values = next_values;
		}

		int queue_updated = 0;
		if (!q.tail ||
		    TRAIL_MAX_SEGMENTS * (nearest_values[0] - q.tail->t) >= trail_fade_t) {
			pos_queue_update(&q, n_fields, nearest_values, trail_fade_t);
			queue_updated = 1;
		}

		if (auto_zoom) {
			const struct pos_it *it = q.head;
#if 0
			double nearest_t, trail_min_t;
			nearest_t = affine_transform(time_scale, nearest_values[0]);
			trail_min_t = affine_inverse_transform(time_scale, nearest_t - 30);
			for (; it && it->t < trail_min_t; it = it->next)
#endif
			center_view3(pt2view, nearest_values+1, n_fields/2,
			             &zctl, display_rate, it, f);
		}

		if (!queue_updated) {
			unsigned i;
			for (i=0; i<n_fields/2; i++) {
				struct pt p = q.tail->pts[i];
				struct pt r = {
					next_values[1+i+i+0],
					next_values[1+i+i+1],
				};
				pt_sub(&r, &p);
				cairo_matrix_transform_distance(pt2view, &r.x, &r.y);
				if (pt_norm_square(&r) >= TRAIL_MAX_DIST * TRAIL_MAX_DIST) {
					pos_queue_update(&q, n_fields, next_values, trail_fade_t);
					break;
				}
			}
		}

		fprintf(stderr, "frame: %lu, rec: %lu, t: %.06e, dmin: (%+5.3e:%+5.3e), dmax: (%+5.3e:%+5.3e) \n",
			frame_nr, rec, cur_values[0], zctl.x0[1], zctl.y0[1], zctl.x1[1], zctl.y1[1]);

		paint_frame(f, pt2view, n_fields, nearest_values, rec, &q,
		            trail_fade_t, (double)rb.valid / rb.sz, 0.0, 1);

		data = cairo_image_surface_get_data(f->sf);
		s = cairo_image_surface_get_stride(f->sf);
		switch (omode) {
		case PPM:
			p = img;
			break;
		case Y4M:
			p = img;
			u = p + f->h * f->w;
			v = u + f->h * f->w / 4;
			break;
		}
		for (y=0; y<f->h; y++, data += s) {
			for (x=0; x<f->w; x++) {
				uint32_t px = ((const uint32_t *)data)[x];
				unsigned r = (px & FRAME_CAIRO_RMASK) >> FRAME_CAIRO_RSHIFT;
				unsigned g = (px & FRAME_CAIRO_GMASK) >> FRAME_CAIRO_GSHIFT;
				unsigned b = (px & FRAME_CAIRO_BMASK) >> FRAME_CAIRO_BSHIFT;
				unsigned char tu, tv;
				switch (omode) {
				case PPM:
					*p++ = r;
					*p++ = g;
					*p++ = b;
					break;
				case Y4M:
					rgb2yuv(p++, &tu, &tv, r, g, b);
#if FRAME_GRAY == 0
#if 1
					unsigned q = (y & 1) << 1 | (x & 1);
					u[x/2] = (q * u[x/2] + tu) / (q+1);
					v[x/2] = (q * v[x/2] + tv) / (q+1);
#else
					if (!(y & 1) && !(x & 1)) {
						u[x/2] = tu / 4;
						v[x/2] = tv / 4;
					} else {
						u[x/2] += tu / 4;
						v[x/2] += tv / 4;
					}
#endif
#endif
					break;
				}
			}
#if FRAME_GRAY == 0
			if (omode == Y4M && (y & 1)) {
				u += f->w / 2;
				v += f->w / 2;
			}
#endif
		}
		switch (omode) {
		case PPM:
			fprintf(fofd, "P6\n%u %u\n255\n", f->w, f->h); /* PPM raw */
			break;
		case Y4M:
			fprintf(fofd, "FRAME\n");
			break;
		}
		fwrite(img, img_sz, 1, fofd);

		frame_nr++;
	}
done:
	free(scratch.buf);
	free(rb.buf);
	free(img);
}

int main(int argc, char **argv)
{
	unsigned n_fields = 0, i;
	const char *tok;
	char *field_spec;
	struct at time_scale = { 0.0, 1.0 };
	struct pt center = { CENTER_X, CENTER_Y };
	double trail_fade_t = -1.0;
	double zoom = SCALE;
	double display_rate = INFINITY;
	int fd = STDIN_FILENO;
	struct frame f;
	int auto_zoom = 0;
	FILE *fofd = stdout;
	enum output_mode omode;

	omode = SDL;

	memset(&f, 0, sizeof(f));
	f.w = WIDTH;
	f.h = HEIGHT;
	f.sdl_flags = SDL_SWSURFACE | SDL_ANYFORMAT | SDL_RESIZABLE;

	static const char *omode_ids[] = {
		[SDL] = "sdl",
		[PPM] = "ppm",
		[Y4M] = "y4m"
	};

	/* option parsing */
	int opt;
	while ((opt = getopt(argc, argv, ":ac:d:f:hi:n:r:s:t:z:")) != -1)
		switch (opt) {
		case 'a': auto_zoom = 1; break;
		case 'c':
			center.x = atof(strtok(optarg, ":"));
			center.y = atof(strtok(NULL, ""));
			break;
		case 'd':
			f.w = atoi(strtok(optarg, "x"));
			f.h = atoi(strtok(NULL, ""));
			break;
		case 'f': trail_fade_t = atof(optarg); break;
		case 'h': goto usage;
		case 'i':
			if (fd != STDIN_FILENO)
				close(fd);
			if ((fd = open(optarg, O_RDONLY)) == -1) {
				perror(optarg);
				exit(1);
			}
			break;
		case 'n':
			if (!strcmp(optarg, omode_ids[SDL])) { omode = SDL; break; }
			if (!strcmp(optarg, omode_ids[PPM])) { omode = PPM; break; }
			if (!strcmp(optarg, omode_ids[Y4M])) { omode = Y4M; break; }
			fprintf(stderr,
				"invalid parameter '%s' for option '-n'\n",
				optarg);
			exit(1);
		case 'r': display_rate = atof(optarg); break;
		case 's': time_scale.s = atof(optarg); break;
		case 't': time_scale.x0 = atof(optarg); break;
		case 'z': zoom = atof(optarg); break;
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
fprintf(stderr, "usage: %s [-OPTS] <t:x0:y0[:x1:y1[:...]]>\n", argv[0]);
fprintf(stderr, "\n");
fprintf(stderr, "Options [defaults in brackets]:\n");
fprintf(stderr, "  -a                      auto-center and -scale view (best used with '-f') [%d]\n", auto_zoom);
fprintf(stderr, "  -c <center-x:center-y>  initial center of view in input units [%g:%g]\n", center.x, center.y);
fprintf(stderr, "  -d <WxH>                initial image dimension in px [%ux%u]\n", f.w, f.h);
fprintf(stderr, "  -f <trail-fade-time>    enable trail of given duration in input time [%g]\n", trail_fade_t);
fprintf(stderr, "  -h                      display this help message\n");
fprintf(stderr, "  -i <input.dat>          use given file as input [stdin]\n");
fprintf(stderr, "  -n { ppm | y4m }        no GUI, write frames as PPM/Y4M-stream on stdout [%s]\n", omode_ids[omode]);
fprintf(stderr, "  -r <display-rate>       update display with max. given Hz / frame-rate [%g]\n", display_rate);
fprintf(stderr, "  -s <time-scale>         scale input timestamps by given factor [%g]\n", time_scale.s);
fprintf(stderr, "  -t <time-offset>        offset input timestamps by given value [%g]\n", time_scale.x0);
fprintf(stderr, "  -z <zoom-factor>        initial scale for view in px/input unit [%g]\n", zoom);
fprintf(stderr, "\n");
fprintf(stderr, "In interactive mode the following actions can be performed:\n");
fprintf(stderr, "  ESC / 'q'               quit the program\n");
fprintf(stderr, "  Arrow l/r/u/d           move view l/r/u/d\n");
fprintf(stderr, "  '+' / '-' / '0'         increase / decrease / reset speed\n");
fprintf(stderr, "  'a'                     toggle auto-center and -scale view\n");
fprintf(stderr, "  Mouse left drag         move view\n");
fprintf(stderr, "  Mouse scroll u/d        zoom in / out\n");
fprintf(stderr, "\n");
fprintf(stderr, "Author: Franz Brauße <mail@franz-brausse.de>, license: GPLv2.\n");
		exit(1);
	}
	field_spec = argv[optind];
	time_scale.x0 *= -time_scale.s;

#if !(defined(_WIN32) || defined(_WIN64))
	/* non-blocking input file to enable interactivity also in case of
	 * waiting for the next input line */
	if (fcntl(fd, F_SETFL, O_NONBLOCK) == -1) {
		fprintf(stderr, "fcntl(%d, F_SETFL, O_NONBLOCK): %s\n",
			fd, strerror(errno));
		exit(1);
	}
#endif

	/* parse field_spec */
	for (tok = field_spec-1; tok; tok = strchr(tok + 1, ':'))
		n_fields++;
	unsigned field_ids[n_fields];
	for (i=0, tok = strtok(field_spec, ":"); tok;
	     i++, tok = strtok(NULL, ":"))
		field_ids[i] = atoi(tok);

	if ((n_fields % 2) != 1) {
		fprintf(stderr, "expected even number of body point samples\n");
		exit(1);
	}

	cairo_matrix_t pt2view;
	cairo_matrix_init_identity(&pt2view);
	cairo_matrix_translate(&pt2view, f.w / 2.0, f.h / 2.0);
	cairo_matrix_scale(&pt2view, zoom, zoom);
	cairo_matrix_translate(&pt2view, -center.x, -center.y);

	switch (omode) {
	case SDL:
		process_interactively(fd, n_fields, field_ids, &f, &pt2view,
		                      &time_scale, trail_fade_t, display_rate,
		                      auto_zoom);
		break;
	case PPM:
	case Y4M:
		process_noninteractively(fd, fofd, omode, n_fields, field_ids,
		                         &f, &pt2view, &time_scale,
		                         trail_fade_t, display_rate, auto_zoom);
		break;
	}
	close(fd);
	fclose(fofd);

	return 0;
}
