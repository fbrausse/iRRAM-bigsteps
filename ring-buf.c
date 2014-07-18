
#include <stdio.h>		/* EOF */
#include <string.h>		/* memchr() */
#include <stdlib.h>		/* malloc(), free() */
#include <assert.h>

#include "ring-buf.h"

void * rb_fill_at(const struct ring_buf *rb, size_t *n)
{
	size_t c1 = rb_c1fill_left(rb);
	size_t c2 = (rb->sz - rb->valid) - c1;
	if (c1) {
		if (n) *n = c1;
		return rb->buf + rb->idx + rb->valid;
	} else {
		if (n) *n = c2;
		return rb->buf + (rb->idx + rb->valid - rb->sz);
	}
}

ssize_t rb_read(struct ring_buf *rb, int fd)
{
	size_t c;
	void *at = rb_fill_at(rb, &c);
	ssize_t r = read(fd, at, c);
	/* fprintf(stderr, "buf read %zd at %zu\n", r, (char *)at - rb->buf); */
	if (r > 0)
		rb->valid += r;
	return r;
}

int rb_getc(struct ring_buf *rb)
{
	if (!rb->valid)
		return EOF;
	char c = rb->buf[rb->idx++];
	if (rb->idx == rb->sz)
		rb->idx = 0;
	rb->valid--;
	return c;
}

int rb_ungetc(struct ring_buf *rb, int c)
{
	if (rb->valid == rb->sz)
		return EOF;
	if (!rb->idx)
		rb->idx = rb->sz;
	rb->buf[--rb->idx] = c;
	rb->valid++;
	return c;
}

int rb_putc(struct ring_buf *rb, int c)
{
	if (rb->valid == rb->sz)
		return EOF;
	*(char *)rb_fill_at(rb, NULL) = c;
	rb->valid++;
	return c & 0xff;
}

void * rb_chr(const struct ring_buf *rb, int c)
{
	size_t c1 = rb_c1drain_left(rb);
	const void *p = rb->buf + rb->idx;
	void *r = memchr(p, c, c1);
	if (r)
		return r;
	return memchr(rb->buf, c, rb->valid - c1);
}

size_t rb_spn(const struct ring_buf *rb, const void *accept, size_t len)
{
	size_t c1, c2;
	const char *p, *p1;
	c1 = rb_c1drain_left(rb);
	p1 = rb->buf + rb->idx;
	c2 = rb->valid - c1;
	for (p = p1; c1; c1--, p++)
		if (!memchr(accept, *p, len))
			return p - p1;
	for (p = rb->buf; c2; c2--, p++)
		if (!memchr(accept, *p, len))
			return rb->sz - (p1 - p);
	return rb->valid;
}

size_t rb_cspn(const struct ring_buf *rb, const void *reject, size_t len)
{
	size_t c1, c2;
	const char *p, *p1;
	c1 = rb_c1drain_left(rb);
	p1 = rb->buf + rb->idx;
	c2 = rb->valid - c1;
	for (p = p1; c1; c1--, p++)
		if (memchr(reject, *p, len))
			return p - p1;
	for (p = rb->buf; c2; c2--, p++)
		if (memchr(reject, *p, len))
			return rb->sz - (p1 - p);
	return rb->valid;
}

int rb_tfer(struct ring_buf *rb2, struct ring_buf *rb1, size_t n)
{
	if (n > rb1->valid)
		return -1;
	if (n > rb2->sz - rb2->valid)
		return -2;
	if (rb1 != rb2)
		assert(rb1->buf != rb2->buf);
	while (n) {
		size_t c1, c2;
		size_t m;
		c1 = rb_c1drain_left(rb1);
		const void *p1 = rb1->buf + rb1->idx;
		void *p2 = rb_fill_at(rb2, &c2);
		m = c1 < c2 ? c1 : c2;
		m = n < m ? n : m;
		memcpy(p2, p1, m);
		n -= m;
		rb_filled(rb2, m);
		rb_drained(rb1, m);
	}
	return 0;
}

int rb_cmp(struct ring_buf rb1, struct ring_buf rb2, size_t n)
{
	while (n) {
		size_t c11 = rb_c1drain_left(&rb1);
		size_t c21 = rb_c1drain_left(&rb2);
		size_t m = c11 < c21 ? c11 : c21;
		int r;
		if (!m)
			return c11 ? -1 : c21 ? 1 : 0;
		m = n < m ? n : m;
		r = memcmp(rb1.buf + rb1.idx, rb2.buf + rb2.idx, m);
		if (r)
			return r;
		n -= m;
		rb_drained(&rb1, m);
		rb_drained(&rb2, m);
	}
	return 0;
}

int rb_memcmp(struct ring_buf rb1, const void *p2, size_t n)
{
	while (1) {
		size_t c1 = rb_c1drain_left(&rb1);
		size_t m = c1 < n ? c1 : n;
		int r;
		if (!m)
			return c1 ? -1 : n ? 1 : 0;
		r = memcmp(rb1.buf + rb1.idx, p2, m);
		if (r)
			return r;
		n -= m;
		rb_drained(&rb1, m);
		p2 = (const char *)p2 + m;
	}
}

int rb_ensure_sz(struct ring_buf *rb, size_t sz)
{
	if (rb->valid >= sz)
		return -2;
	if (rb->sz >= sz)
		return 0;
	struct ring_buf tmp = { malloc(sz), 0, 0, sz };
	if (!tmp.buf)
		return -1;
	int ret = rb_tfer(&tmp, rb, rb->valid);
	(void)ret;
	assert(ret == 0);
	free(rb->buf);
	*rb = tmp;
	return 0;
}
