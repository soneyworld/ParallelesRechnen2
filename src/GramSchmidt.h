/*
 * GramSchmidt.h
 *
 *      Author: Till Lorentzen
 */

#ifndef GRAMSCHMIDT_H_
#define GRAMSCHMIDT_H_
#define NUMBER double

NUMBER dot(NUMBER *a, NUMBER *b, unsigned int length);

void gs(NUMBER *v,NUMBER *w, unsigned int length, unsigned int lda);

#endif /* GRAMSCHMIDT_H_ */
