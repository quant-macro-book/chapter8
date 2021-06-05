## 『経済セミナー』「定量的マクロ経済学と数値計算」

#### 連載第8回(2・3月号)

* **2 (クルセル=スミスモデル)** の結果を再現するファイル -> 内挿法(Interpolation)に関してはNumerical Recipesのサブルーチンを利用しています。そのため、ダウンロードしただけではコンパイルは通りません。下記のサブルーチンについて、各自で補完してください。
  * locate, polint, splint, spline, tridag_ser, assert_eq2, assert_eq3, assert_eq4, assert_eqn, iminloc, nrerror
  * Matlabファイルについてもlocate.mが含まれていないため、途中でエラーが出ます。同じく、Numerical Recipesのlocate.f90をMatlabに翻訳したものを各自で追加してください。
  * オイラー方程式を解く際にEndogenous Gridpoint Method(EGM)という方法を使っています。EGMに関する簡単なメモは、MATLAB_KSの中に追加してあります。
  * Press et al. (1996) ``Numerical Recipes in Fortran 90," Cambridge University Press.
* **3.2 (HANKモデル)定常状態の計算**の結果を再現するファイル -> main_ti_cheb.m, nti_cheb.m, EulerEq_cheb.m, calcerr.m, CRRA.m, mu_CRRA.m, polygrid.m, polybas.m

#### 注意
* 2節(クルセル=スミスモデル)のMATLABコード、Fortranコードについては山田知明が作成しています。
* 3節(HANKモデル)のMATLABコード、Juliaコードについては砂川武貴が作成しています。
<!--* 文字コードがUTF-8のため、一部の日本語がWindowsでは正しく表示されない可能性があります(to be fixed.)
#### 未完成
* Pythonコードはこれからアップ予定です。
