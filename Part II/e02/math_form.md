# Mathematical Reformulation

$$
G(z) = \frac{\theta_3z^2 + \theta_4z + \theta_5}{z^2 + \theta_1z + \theta_2},
\;\;\;\;
G(z) = \frac{Y(z)}{U(z)}
$$
$$
Y(z) = \frac{\theta_3z^2 + \theta_4z + \theta_5}{z^2 + \theta_1z + \theta_2} \, U(z)
$$
$$
Y(z) = \frac{\theta_3 + \theta_4z^{-1} + \theta_5z^{-2}}{1 + \theta_1z^{-1} + \theta_2z^{-2}} \, U(z)
$$
$$
\left( 1 + \theta_1z^{-1} + \theta_2z^{-2} \right) Y(z) = \left( \theta_3 + \theta_4z^{-1} + \theta_5z^{-2} \right) U(z)
$$
$$
Y(z) + \theta_1z^{-1}Y(z) + \theta_2z^{-2}Y(z) - \theta_3U(z) - \theta_4z^{-1}U(z) - \theta_5z^{-2}U(z) = 0
$$
$$
y(k) + \theta_1q^{-1}y(k) + \theta_2q^{-2}y(k) - \theta_3u(k) - \theta_4q^{-1}u(k) - \theta_5q^{-2}u(k) = 0
$$
$$
y(k) + \theta_1y(k-1) + \theta_2y(k-2) - \theta_3u(k) - \theta_4u(k-1) - \theta_5u(k-2) = 0
$$
$$
\tilde{y}(k) - \eta(k) + \theta_1\tilde{y}(k-1) - \theta_1\eta(k-1) + \theta_2\tilde{y}(k-2) - \theta_2\eta(k-2) + \\
- \theta_3u(k) - \theta_4u(k-1) - \theta_5u(k-2) = 0
$$
$$
\mathcal{D}_{\theta} = \{ \theta \in \mathbb{R}^p : \\
\tilde{y}(k) - \eta(k) + \theta_1\tilde{y}(k-1) - \theta_1\eta(k-1) + \theta_2\tilde{y}(k-2) - \theta_2\eta(k-2) + \\
- \theta_3u(k) - \theta_4u(k-1) - \theta_5u(k-2) = 0, \;\; \forall k = 3, \dots, N,\\
\eta(k) \in \mathcal{B}_{\eta} = \{ \left| \eta(k) \right| \le \Delta\eta(k) \} \}
$$
$$
\mathcal{D}_{\theta, \eta} = \{ \theta \in \mathbb{R}^{p=5}, \eta \in \mathbb{R}^{N=50} : \\
\tilde{y}(k) - \eta(k) + \theta_1\tilde{y}(k-1) - \theta_1\eta(k-1) + \theta_2\tilde{y}(k-2) - \theta_2\eta(k-2) + \\
- \theta_3u(k) - \theta_4u(k-1) - \theta_5u(k-2) = 0, \;\; \forall k = 3, \dots, N, \\
-5 = -\Delta\eta \le \eta \le \Delta\eta = 5, \;\; \forall k = 1,\dots, N\}
$$
$$
PUI_{\theta_i} = \left[ \underline{\theta}_i, \bar{\theta}_i \right], \;\;\;\; i = 1,\dots,5
$$
$$
\underline{\theta}_i = \min\limits_{\theta,\eta \in \mathcal{D}_{\theta,\eta}} \theta_i, \;\;\;\; \bar{\theta}_i = \min\limits_{\theta,\eta \in \mathcal{D}_{\theta,\eta}} \left( -\theta_i \right)
$$
$$
\theta_c = \min\limits_{\theta \in \mathbb{R}^5} \max\limits_{\theta' \in \mathcal{D}_{\theta}} \left|\left| \theta - \theta' \right|\right|_{\infty}
$$
$$
\theta_{c,i} = \frac{\bar{\theta}_i + \underline{\theta}_i}{2}
$$

- For translating into code
$$
1) \;\; \tilde{y}(k)    \\
2) \;\; - \eta(k)   \\
3) \;\; \theta_1\tilde{y}(k-1)  \\
4) \;\; - \theta_1\eta(k-1) \\
5) \;\; \theta_2\tilde{y}(k-2)  \\
6) \;\; - \theta_2\eta(k-2) \\
7) \;\; - \theta_3u(k)  \\
8) \;\; - \theta_4u(k-1)    \\
9) \;\; - \theta_5u(k-2)
$$

- To force the dc-gain to be equal to a certain value (e.g. $118$), one should add another constraint evaluating the transfer function for $z \rightarrow 1$
$$
G(z) = \frac{\theta_3z^2 + \theta_4z + \theta_5}{z^2 + \theta_1z + \theta_2}
$$
$$
\lim\limits_{z \rightarrow 1} G(z) = \lim\limits_{z \rightarrow 1} \frac{\theta_3z^2 + \theta_4z + \theta_5}{z^2 + \theta_1z + \theta_2} =
\frac{\theta_3 + \theta_4 + \theta_5}{1 + \theta_1 + \theta_2} = 118
$$
$$
\theta_3 + \theta_4 + \theta_5 = 118 \left( 1 + \theta_1 + \theta_2 \right)
$$
$$
118 + 118\theta_1 + 118\theta_2 - \theta_3 - \theta_4 - \theta_5 = 0
$$