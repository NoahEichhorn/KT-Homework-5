# Homework 5 - Physical parameters of a Dirac spinors 

by Hanna Schulteis (108023203549), Noah Eichhorn (108020260079), Nele Blume (108023203682)


## Spinor für Nele Blume

   Phi_1 = 14 e^(2i)  
   Phi_2 = 5 e^(12i)  
   Phi_3 = 12 e^(21i)  
   Phi_4 = 5 e^(13i)  

Dieser Spinor ist ungültig, da die Relation E^2 = abs(vec(p))^2 + m^2 nicht gilt.  
Mit den auf dem Arbeitsblatt angegebenen Formeln für E und m, und mit 2p_i = Phi^dagger gamma^0 gamma^i Phi lässt sich E=195, m=26 und abs(vec(p))^2=35213,07 berechnen. 
Um den Spinor durch Änderung von höchstens zwei Komponenten in einen gültigen Spinor umzuändern, tauschen wir die letzten beiden Komponenten aus:

   Phi_1^Strich = 14 e^(2i)  
   Phi_2^Strich =  5 e^(12i)  
   Phi_3^Strich = 14 e^(2i)  
   Phi_4^Strich = -5 e^(12i)  

Es wird für diesen nun gültigen Spinor berechnet: E = 221, m = 0, vec(p) = (0,0,221), vec(S)_0 = (-0,2658; -0,1723; 0,3869)  
Zur besseren Lesbarkeit können die Rechnungen (für alle Namen) dazu im Dokument `Rechnungen_HA5.pdf` eingesehen werden.

Als nächstes suchen wir den gültigen Spinor mit dem kleinsten Abstand zum Originalen. Dafür definieren wir den Impuls in Kugelkoordinaten:

   vec(p) = p(sin(theta)cos(phi), sin(theta)sin(phi), cos(theta))

Außerdem parametrisieren wir das Chi, also den 2-Spinor aus der Form, in die wir den Spinor bringen wollen, wie folgt: 

   Chi(alpha, beta) = (cos(alpha/2), e^(i beta)sin(alspha/2))
   
So erreichen wir, dass E^2 = p^2 + m^2 für jeden Satz an Parametern gilt. Der Abstand hat also die Abhängigkeit d(m,p,theta,phi,alpha,beta). Da dies ein Minimierungsproblem mit 5 Parametern ist, haben wir uns für eine numerische Lösung entschieden.  
Speziell für den Namen "Nele Blume" macht das der Code `closest spinor nele.py`. Dieser Code macht einmal die Berechnungen für einen u und einmal für einen v Spinor und entscheidet dann, welcher näher ist.

Der Output ist: 

Originaler Spinor φ:  
  φ[1] = -5.826056 + 12.730164 i  
  φ[2] = 4.219270 - 2.682865 i  
  φ[3] = -6.572751 + 10.039868 i  
  φ[4] = 4.537234 + 2.100835 i  

||φ||^2 = 390.0  

Starte numerische Minimierung für u-Spinoren (BFGS) ...  
u-Minimierung erfolgreich? False  
Minimale Distanz d_u_min = 2.7779488934502865  

Starte numerische Minimierung für v-Spinoren (BFGS) ...  
v-Minimierung erfolgreich? False  
Minimale Distanz d_v_min = 3.674273433223334  

===== Optimaler u-Spinor =====  
  u_opt[1] = -7.136134 + 12.115104 i  
  u_opt[2] = 4.024885 - 2.168762 i  
  u_opt[3] = -5.116323 + 10.897734 i  
  u_opt[4] = 4.276512 + 1.616980 i  

Physikalische Größen (u):  
  E_u = 192.2219665221858  
  m_u = 26.381421812714844  
  p_u = [-55.1550771  -96.11540724 154.83233176]  
  ||u_opt||^2 = 384.44393304437335  
  Distanz d_u (Definition)         = 2.777948893450372  
  Distanz d_u (Tabellenkonvention) = 5.555897786900744  

Ruhe-Spin S^0 (u):  
  Sx^0 = -0.25158298754552916  
  Sy^0 = -0.15226352964142958  
  Sz^0 = 0.4043783103960938  
  |S^0| = 0.5  

===== Optimaler v-Spinor =====  
  v_opt[1] = -6.772860 + 11.288112 i  
  v_opt[2] = 3.750872 - 1.990803 i  
  v_opt[3] = -5.589152 + 11.649419 i  
  v_opt[4] = 4.610023 + 1.767925 i  

Physikalische Größen (v):  
  E_v = 191.32545227261076  
  m_v = 2.166477025263399e-09  
  p_v = [-55.42239658 -96.58095447 155.58247292]  
  ||v_opt||^2 = 382.6509045452235  
  Distanz d_v (Definition)         = 3.674273433223338  
  Distanz d_v (Tabellenkonvention) = 7.348546866446676  

Ruhe-Spin S^0 (v):  
  Sx^0 = 0.027026313930272456  
  Sy^0 = 0.3323410197457544  
  Sz^0 = -0.37258425214935537  
  |S^0| = 0.5  

===== Vergleich =====  
  d_u (Def.) = 2.777948893450372  
  d_v (Def.) = 3.674273433223338  

Ergebnis: φ ist näher an einem u-Spinor als an einem v-Spinor.  

Process finished with exit code 0


Der nächste Spin ist also der angegebene optimale u-Spinor mit den angegebenen physikalischen Größen.




## Spinor für Hanna Schulteis

   Phi_1 = 8 e^(19i)  
   Phi_2 = 1 e^(3i)  
   Phi_3 = 14 e^(8i)  
   Phi_4 = 14 e^(21i)  

Auch dieser Spinor ist ungültig, weil man(, wie bei dem Fall für Nele Blume,) E^2 = 52212,25; m^2 = 26732,25; abs(vec(p))^2 = 15180,78 erhält. Damit ist E^2 = m^2 + abs(vec(p))^2 nicht erfüllt.

Wir verändern ihn mit 2 Komponenten, sodass er gültig ist:

   Phi_1^Strich = 8 e^(19i)  
   Phi_2^Strich =  1 e^(3i)  
   Phi_3^Strich = 8 e^(19i)  
   Phi_4^Strich = -1 e^(3i)  

Dafür erhalten wir:  
   E = 65; m = 0; ; vec(p) = (0, 0, 65); abs(vec(p)) = 65.
   
Es ist also E^2 = m^2 + abs(vec(p))^2 = 4225 erfüllt.  
Der Ruhespin ist nach demselben Prinzip:

   vec(S)^0 = (-0,118; 0,035; 0,485)

Für den allernächsten Spin erhält man den folgenden Output: 

Erzeugter Spinor φ aus Namen:  
  φ[1] = 7.909637 + 1.199018 i  
  φ[2] = -0.989992 + 0.141120 i  
  φ[3] = -2.037000 + 13.851015 i  
  φ[4] = -7.668210 + 11.713179 i  

||φ||^2 = 457.0  
x_vals (Vorname-Positionen) = [8, 1, 14, 14]  
y_vals (Nachname-Positionen) = [19, 3, 8, 21]  

Starte numerische Minimierung für u-Spinoren (BFGS) ...  
u-Minimierung erfolgreich? False  
Minimale Distanz d_u_min = 52.64484965552009  

Starte numerische Minimierung für v-Spinoren (BFGS) ...  
v-Minimierung erfolgreich? False  
Minimale Distanz d_v_min = 11.886642058080383  

===== Optimaler u-Spinor =====  
  u_opt[1] = 10.832909 + 1.667969 i  
  u_opt[2] = -6.893765 - 2.862966 i  
  u_opt[3] = -1.062015 + 7.321580 i  
  u_opt[4] = -5.798672 + 9.353962 i  

Physikalische Größen (u):  
  E_u = 175.85461946219917  
  m_u = 1.6534847490814774e-09  
  p_u = [-60.85451599 164.5164015  -12.48714169]  
  ||u_opt||^2 = 351.7092389244004  
  Distanz d_u (Definition)       = 52.644849655520105  
  Distanz d_u (Tabellenkonvention) = 105.28969931104021  

Ruhe-Spin S^0 (u):  
  Sx^0 = -0.4518213630988246  
  Sy^0 = -0.1109761202826331  
  Sz^0 = 0.1831440869221135  
  |S^0| = 0.5  

===== Optimaler v-Spinor =====  
  v_opt[1] = 4.935831 + 0.763753 i  
  v_opt[2] = -3.950639 - 1.751290 i  
  v_opt[3] = -1.858089 + 12.694784 i  
  v_opt[4] = -8.078460 + 12.638638 i  

Physikalische Größen (v):  
  E_v = 216.6135141456115  
  m_v = 172.99320554504882  
  p_v = [-45.11271989 121.95868666  -9.25670558]  
  ||v_opt||^2 = 433.2270282912242  
  Distanz d_v (Definition)       = 11.886642058080394  
  Distanz d_v (Tabellenkonvention) = 23.773284116160788  

Ruhe-Spin S^0 (v):  
  Sx^0 = -0.45033945672504716  
  Sy^0 = -0.2029497498455389  
  Sz^0 = 0.07749692093381831  
  |S^0| = 0.5  

===== Vergleich =====  
  d_u (Def.) = 52.644849655520105  
  d_v (Def.) = 11.886642058080394  

Ergebnis: φ ist näher an einem v-Spinor als an einem u-Spinor.

Process finished with exit code 0




## Spinor für Noah Eichhorn

   Phi_1 = 14 e^(5i)  
   Phi_2 = 15 e^(9i)  
   Phi_3 = 1 e^(3i)  
   Phi_4 = 8 e^(8i)  

Für diesen Spinor ist E^2 = 59049; m^2 = 31684, abs(vec(p))^2 = 14435,79. Auch hier ist also E^2 = m^2 + p^2 nicht erfüllt.


Wir verändern ihn mit 2 Komponenten, sodass er gültig ist:

   Phi_1^Strich = 14 e^(5i)  
   Phi_2^Strich =  15 e^(9i)  
   Phi_3^Strich = 14 e^(5i)  
   Phi_4^Strich = -15 e^(9i)  

Es ist E = 421; m = 0; vec(p) = (0, 0, 421); abs(vec(p)) = 421; vec(S)^0 = (-0,326; -0,378, -0,034)  
E^2 = m^2 + abs(vec(p))^2 = 421^2 ist also erfüllt.

Für den nächsten Spinor bekommen wir den Output: 


Erzeugter Spinor φ aus Namen:  
  φ[1] = 3.971271 - 13.424940 i  
  φ[2] = -13.666954 + 6.181777 i  
  φ[3] = -0.989992 + 0.141120 i  
  φ[4] = -1.164000 + 7.914866 i  

||φ||^2 = 486.0  
x_vals (Vorname-Positionen) = [14, 15, 1, 8]  
y_vals (Nachname-Positionen) = [5, 9, 3, 8]  

Starte numerische Minimierung für u-Spinoren (BFGS) ...  
u-Minimierung erfolgreich? False  
Minimale Distanz d_u_min = 14.122406394779603  

Starte numerische Minimierung für v-Spinoren (BFGS) ...  
v-Minimierung erfolgreich? False  
Minimale Distanz d_v_min = 61.42548287258836  

===== Optimaler u-Spinor =====  
  u_opt[1] = 4.269780 - 14.045681 i  
  u_opt[2] = -12.470357 + 6.896447 i  
  u_opt[3] = 2.498957 + 1.201737 i  
  u_opt[4] = -2.877082 + 4.817211 i  

Physikalische Größen (u):  
  E_u = 228.8775193321626  
  m_u = 189.705438616785  
  p_u = [-102.82068961   12.37786952  -75.30909341]  
  ||u_opt||^2 = 457.7550386643264
  Distanz d_u (Definition)       = 14.122406394779677  
  Distanz d_u (Tabellenkonvention) = 28.244812789559354  

Ruhe-Spin S^0 (u):  
  Sx^0 = -0.35861702344282287  
  Sy^0 = -0.3480990826862359  
  Sz^0 = 0.014861330021601152  
  |S^0| = 0.5  

===== Optimaler v-Spinor =====  
  v_opt[1] = 3.126615 - 9.875426 i  
  v_opt[2] = -6.785110 + 5.313819 i  
  v_opt[3] = 4.123099 + 2.196956 i  
  v_opt[4] = -5.546471 + 11.357124 i  

Physikalische Größen (v):  
  E_v = 181.57415309044617  
  m_v = 5.418104310247578e-09  
  p_v = [-145.79956239   17.55157701 -106.78765295]  
  ||v_opt||^2 = 363.14830618089434  
  Distanz d_v (Definition)       = 61.42548287258833  
  Distanz d_v (Tabellenkonvention) = 122.85096574517667  

Ruhe-Spin S^0 (v):  
  Sx^0 = -0.011468902206465275  
  Sy^0 = -0.3250016296568103  
  Sz^0 = 0.37979258155287354  
  |S^0| = 0.5  

===== Vergleich =====  
  d_u (Def.) = 14.122406394779677  
  d_v (Def.) = 61.42548287258833  

Ergebnis: φ ist näher an einem u-Spinor als an einem v-Spinor.

Process finished with exit code 0




## Bonusaufgabe

Die Datei `closest spinor name.py` nimmt Vor- und Nachnamen mit jeweils mind. vier Buchstaben (für kürzere Namen war keine Angabe gegeben, wie der Spinor aussehen soll) entgegen und gibt den nächsten u- und v-Spinor mit den jeweiligen physikalischen Größen aus. Dazu entscheidet es, welcher von beiden näher am ursprünglichen Spinor liegt.

Für den Cross-Check erhält man für “Mikhail Mikhasenko”:

Erzeugter Spinor φ aus Namen:  
  φ[1] = 11.796808 + 5.462171 i  
  φ[2] = -8.200172 + 3.709066 i  
  φ[3] = 0.048683 - 10.999892 i  
  φ[4] = -1.164000 + 7.914866 i  

||φ||^2 = 435.0  
x_vals (Vorname-Positionen) = [13, 9, 11, 8]  
y_vals (Nachname-Positionen) = [13, 9, 11, 8]  

Starte numerische Minimierung für u-Spinoren (BFGS) ...  
u-Minimierung erfolgreich? False  
Minimale Distanz d_u_min = 56.37633554980869  

Starte numerische Minimierung für v-Spinoren (BFGS) ...  
v-Minimierung erfolgreich? False  
Minimale Distanz d_v_min = 58.96106461101567  

===== Optimaler u-Spinor =====  
  u_opt[1] = 8.137397 + 8.357696 i  
  u_opt[2] = -5.412016 + 6.764445 i  
  u_opt[3] = -4.895057 - 6.186488 i  
  u_opt[4] = -5.165304 + 4.713569 i  

Physikalische Größen (u):  
  E_u = 161.12412796941314  
  m_u = 49.991809789696816  
  p_u = [ -17.99358484   14.93249393 -151.37719477]  
  ||u_opt||^2 = 322.2482559388278  
  Distanz d_u (Definition)       = 56.37633554980869  
  Distanz d_u (Tabellenkonvention) = 112.75267109961737  

Ruhe-Spin S^0 (u):  
  Sx^0 = 0.05918763659694054  
  Sy^0 = 0.4749852563196924  
  Sz^0 = 0.1445193065060327  
  |S^0| = 0.5  

===== Optimaler v-Spinor =====  
  v_opt[1] = 6.328463 + 7.757035 i  
  v_opt[2] = -4.142154 + 6.415244 i  
  v_opt[3] = -5.141220 - 8.017058 i  
  v_opt[4] = -5.594515 + 6.043994 i  

Physikalische Größen (v):  
  E_v = 158.53383417678128  
  m_v = 1.4141203152039095e-09  
  p_v = [ -18.65686876   15.45597186 -156.67166547]  
  ||v_opt||^2 = 317.06766835356444  
  Distanz d_v (Definition)       = 58.96106461101571  
  Distanz d_v (Tabellenkonvention) = 117.92212922203142  

Ruhe-Spin S^0 (v):  
  Sx^0 = 0.12421586860610012  
  Sy^0 = 0.4789202238073532  
  Sz^0 = -0.0721514879593402  
  |S^0| = 0.5  

===== Vergleich =====  
  d_u (Def.) = 56.37633554980869  
  d_v (Def.) = 58.96106461101571  

Ergebnis: φ ist näher an einem u-Spinor als an einem v-Spinor.

Process finished with exit code 0

Es fällt auf, dass der Abstand mit der auf dem Blatt gegebenen Defintion genau die Hälfte desjenigen im Cross-Check liefert, daher wurde der Abstand in der Tabellenkonvention eingefügt.  
Für  “Misha Mikhasenko” erhält man analog:
 
 Erzeugter Spinor φ aus Namen:  
  φ[1] = 11.796808 + 5.462171 i  
  φ[2] = -8.200172 + 3.709066 i  
  φ[3] = 0.084088 - 18.999814 i  
  φ[4] = -1.164000 + 7.914866 i  

||φ||^2 = 675.0  
x_vals (Vorname-Positionen) = [13, 9, 19, 8]  
y_vals (Nachname-Positionen) = [13, 9, 11, 8]  

Starte numerische Minimierung für u-Spinoren (BFGS) ...  
u-Minimierung erfolgreich? False  
Minimale Distanz d_u_min = 89.81804702924819  

Starte numerische Minimierung für v-Spinoren (BFGS) ...  
v-Minimierung erfolgreich? False  
Minimale Distanz d_v_min = 78.50413749301026  

===== Optimaler u-Spinor =====  
  u_opt[1] = 4.616381 + 10.007759 i  
  u_opt[2] = -7.989173 + 7.898626 i  
  u_opt[3] = -4.825123 - 13.888837 i  
  u_opt[4] = -4.853875 + 2.817745 i  

Physikalische Größen (u):  
  E_u = 247.68139456360834  
  m_u = 2.3051146996811156e-10  
  p_u = [ -65.36201211  -87.48796411 -222.30550313]  
  ||u_opt||^2 = 495.36278912721866  
  Distanz d_u (Definition)       = 89.8180470292483  
  Distanz d_u (Tabellenkonvention) = 179.6360940584966  

Ruhe-Spin S^0 (u):  
  Sx^0 = 0.1702448356073134  
  Sy^0 = 0.47002636963148925  
  Sz^0 = -0.009586855588850758  
  |S^0| = 0.5  

===== Optimaler v-Spinor =====  
  v_opt[1] = 1.917615 + 7.771463 i  
  v_opt[2] = -5.513908 + 6.241861 i  
  v_opt[3] = -4.194525 - 17.944104 i  
  v_opt[4] = -4.600483 + 4.879170 i  

Physikalische Größen (v):  
  E_v = 258.9962733829825  
  m_v = 125.55937808318313  
  p_v = [ -59.78005332  -80.01492383 -203.31716453]  
  ||v_opt||^2 = 517.9925467659664  
  Distanz d_v (Definition)       = 78.50413749301038  
  Distanz d_v (Tabellenkonvention) = 157.00827498602075  

Ruhe-Spin S^0 (v):  
  Sx^0 = 0.1774918528148806  
  Sy^0 = 0.2678867050068684  
  Sz^0 = -0.383057900930008  
  |S^0| = 0.5  

===== Vergleich =====  
  d_u (Def.) = 89.8180470292483  
  d_v (Def.) = 78.50413749301038  

Ergebnis: φ ist näher an einem v-Spinor als an einem u-Spinor.

Process finished with exit code 0
