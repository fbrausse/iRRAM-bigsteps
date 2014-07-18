
#ifndef RING_BUF_H
#define RING_BUF_H

#include <unistd.h>				/* size_t */

#if __STDC_VERSION__ >= 199901L
# define INLINE			inline
#else
# define INLINE
#endif

struct ring_buf {
	char *buf;
	size_t idx, valid, sz;
};

/* Generally this buffer is in one of the two states:
 *
 *   lower          buf                          buf
 * addresses      +------+                     +------+
 *                |      |                     | #### |
 *     |          +------+ <-- idx             +------+ <-- (idx+valid)-sz
 *     |          | #### |                     |      |
 *     |          | #### |                     |      |
 *     v          | #### |                     |      |
 *                +------+ <-- idx+valid       +------+ <-- idx
 *   higher       |      |                     | #### |
 * addresses      +------+ <-- sz              +------+ <-- sz
 */

#define RING_BUF_INIT		{ NULL, 0, 0, 0 }
#define RING_BUF_WRAP(ptr,sz)	(struct ring_buf){ (void *)(ptr),0,(sz),(sz) }

/* Returns the number of valid bytes in chunk 1, that is beginning at rb->idx.
 * This value always is <= rb->valid; if rb->valid != 0 it is >= 1. */
static INLINE size_t rb_c1drain_left(const struct ring_buf *rb)
{
	if (rb->idx + rb->valid >= rb->sz)
		return rb->sz - rb->idx;
	return rb->valid;
}

static INLINE size_t rb_c1fill_left(const struct ring_buf *rb)
{
	if (rb->idx + rb->valid >= rb->sz)
		return 0;
	return rb->sz - (rb->idx + rb->valid);
}

static INLINE int rb_drained(struct ring_buf *rb, size_t n)
{
	if (n > rb->valid)
		return -1;
	rb->valid -= n;
	rb->idx   += n;
	if (rb->idx >= rb->sz)
		rb->idx -= rb->sz;
	return 0;
}

static INLINE int rb_filled(struct ring_buf *rb, size_t n)
{
	if (rb->valid + n > rb->sz)
		return -1;
	rb->valid += n;
	return 0;
}

static INLINE void rb_reset(struct ring_buf *rb)
{
	rb->idx = rb->valid = 0;
}

static INLINE size_t rb_ptrdiff(const struct ring_buf *rb, const void *p)
{
	size_t q = (const char *)p - rb->buf;
	return (q < rb->idx ? rb->sz : 0) + q - rb->idx;
}

void *  rb_fill_at(const struct ring_buf *rb, size_t *n);
ssize_t rb_read(struct ring_buf *rb, int fd);
int     rb_getc(struct ring_buf *rb);
int     rb_ungetc(struct ring_buf *rb, int c);
int     rb_putc(struct ring_buf *rb, int c);
void *  rb_chr(const struct ring_buf *rb, int c);
size_t  rb_spn(const struct ring_buf *rb, const void *accept, size_t len);
size_t  rb_cspn(const struct ring_buf *rb, const void *reject, size_t len);
int     rb_tfer(struct ring_buf *rb2, struct ring_buf *rb1, size_t n);
int     rb_ensure_sz(struct ring_buf *rb, size_t sz);
int     rb_cmp(struct ring_buf rb1, struct ring_buf rb2, size_t n);
int     rb_memcmp(struct ring_buf rb1, const void *p2, size_t n);

#endif
