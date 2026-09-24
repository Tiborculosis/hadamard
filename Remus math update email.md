Hi Remus,

Organizing my thoughts has never been a strength of mine. But I'll try to put my main results together as cleanly as possible. Much of this may be redundant since you're aware of my earlier results (and responsible for some of the preliminaries).

The first building blocks of these results would be the uncertainty principles. All of these somehow relate bounding the support of a vector $x$ in terms of its relationship to some matrix $U$ with certain properties. For convenience, let's say that an $n\times n$ matrix $U$ is $k$-*fancy* if it has no $k\times k$ submatrix whose kernel contains a vector $v$ which is supported on $k$ entries, and that $U$ is *fancy* if it is $k$-fancy for all $k\leq n$. Also, let $D_x$ represent $\operatorname{diag}(x)$ in all these statements.

# Uncertainty principles
### $x$ and $Ux$
* Let $U\in\mathbb{C}^{n\times n}$.
  For every nonzero $x\in\mathbb{C}^n$, define $S=\operatorname{supp}(x)$ and $T=\operatorname{supp}(Ux)$ with $s=|S|$ and $t=|T|$.
  Then, $x_S\in\ker(U_{J\times S}\Leftrightarrow J\subseteq T^c$.
  Moreover, if $n-t\geq s$, then $x_S$ is in the kernel of an $s\times s$ minor of $U$.
### $x$ and the columns of $UD_xU^*$
* Let $U\in\mathbb{C}^{n\times n}$. For every nonzero $x\in\mathbb{C}^n$, define $\widehat{x}_k$ as the $k^{\text{th}}$ column of $UD_xU^*$.
  Next, define $S=\operatorname{supp}(x)$ and $T_k=\operatorname{supp}(\widehat{x}_k)$ , with $s=|S|$ and $t_k=|T_k|$.
  Then, $(D_x)_{S\times S}U^*_{S\times k}\in\ker(U_{J\times S})\Leftrightarrow J\subseteq T_k^c$.
  Moreover, if $n-t_k\geq s$, then $U$ is not $s$-fancy.
### $x$ and the columns of $UD_xU^*a$
* Let $U,a\in\mathbb{C}^{n\times n}$, and let $x\in\mathbb{C}^n$ be nonzero with $\widehat{x}_k^a$ the $k^{\text{th}}$ column of $(UD_xU^*a)$.
  Next, define $S=\operatorname{supp}(x)$ with $s=|S|$.
  Consider the matrix $Z=(D_x)_{S\times S}(U^*a)_{S\times\{1,\ldots, n\}}$, with $z_k$ the $k^{\text{th}}$ column of $Z$.
  Further define $R_k=\operatorname{supp}(z_k)$ and $T_k=\operatorname{supp}(\widehat{x}_k^a)$ with $r_k=|R_k|$ and $t_k=|T_k|$.
  Then, for every $J\subseteq\{1,\ldots,n\}$, $U_{J\times R_k}(z_k)_{R_k}=0\Leftrightarrow J\subseteq T_k^c$.
  Moreover, if $z_k\neq0$ and $n-t_k\geq r_k$, then $U$ is not $r_k$-fancy.
  If $U$ and $a$ are invertible, then $Z$ has rank $s$, so $z_k\neq0$ for at least $s$ distinct values of $k$.
  Hence, at least $s$ columns of $UD_xU^*a$ have support of size at least $n-s+1$.

The proofs of the uncertainty principles are all fairly accessible, but I can send them if needed. Importantly, if $U$ is fancy, then all of these uncertainty principles are strong enough to be used in a Haagerup-like argument.

Speaking of Haagerup-like arguments,
# Haagerup-like arguments
### $UDU^*$ Hadamards
* If $U$ is fancy, then there are at most finitely-many Hadamard matrices of the form $UDU^*$.
* This is one of your results, but extended from $U$ superregular to $U$ fancy.
### $UDU^*$-core Hadamards
* If $U$ is fancy, then there are at most finitely-many Hadamard matrices that are $UDU^*$-core.
* This is a slight generalization of the result I presented at GPOTS in May.
### $UDU^*$ $a$-Hadamards that commute with $a$
* If $U$ is fancy, then at most finitely many matrices of the form $UDU^*$ $a$-Hadamard for some $a$ with which they commute.
* The argument generalizes nicely from the Hadamard case with the new uncertainty principle.
### $UDU^*$-core $a$-Hadamards that commute with $a$
* If $U$ is fancy,  then there are at most finitely many $UDU^*$-core matrices that are $a$-Hadamard for some $a$ with which they commute.
* This argument generalizes nicely from the Hadamard case with the new uncertainty principle.
# $UDU^*$-core generalization

Given a vector $x$, let $M(x)$ to be the $C(x)$-core $n+1\times n+1$ matrix, where $C(x)=UD_xU^*$.
Taking $U$ and $M(x)$ unitary is enough to prove $x_1=-1$ and $|x_k|=\sqrt{n+1}$ for all $k\in\{2,\ldots, n\}$. Furthermore, we may conclude that all $UDU^*$-core *unitaries* commute with one another.
We can actually explicitly write the unitary that generates that MASA for any given $n$.
The variety conditions $x_1=y_1=1$ and $x_ky_k=n+1$ for $k\in\{2,\ldots, n\}$ hold from these alone.
$M(x)$ only needs to be Hadamard to get $C(x)_{j,k}C(y)_{k,j}=1$ for all $j,k\in\{1,\ldots, n\}$.
And $U$ only needs to be fancy to draw the final contradiction with the uncertainty principle.

There's a good chance I'm forgetting something important, but this seems like enough results to warrant sending, even if I think of something I'm forgetting later. I haven't had much opportunity since Thursday to explore the other notions of near-Hadamard matrices, but I have actually been finding a decent amount of momentum in the Haagerup-like generalizations. If you'd like me to put those on hold while I explore the other notions of near-Hadamard matrices then I can, but I don't know how long it might take me to find the same sort of inspiration.

Alternatively, I can try to explore the generalizations we talked about in your office last week. Specifically regarding the $a$-Hadamard case, trying to figure out what the *actual* condition needs to be, since regular Hadamard matrices are also $a$-Hadamard for diagonal $a$, but don't necessarily commute with said matrices. I'm not quite sure how to progress in this direction, but it *does* seem like a natural point of inquiry. Let me know what you think and I'll do my best to continue making progress towards a dissertation.

Best,
Tibor