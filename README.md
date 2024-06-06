## 北尾早霧・砂川武貴・山田知明『定量的マクロ経済学と数値計算』日本評論社

### 第8章：定量的マクロ経済学のフロンティアを覗いてみる

#### Julia
* 各節の番号に対応したフォルダにコードが格納されています。
* `8_2_Krusell_Smith`：RAの小野泰輝さんに作成していただいたJulia用再現コードです。もちろん残りうるあらゆる間違いは、すべて著者たちの責任です。
* `8_4_HANK`：8.4節の結果を再現するためのコード。

#### MATLAB
* 各節の番号に対応したフォルダにコードが格納されています。
* `8_2_Krusell_Smith`：`locate.m`が含まれていないため、途中でエラーが出ます。[Numerical Recipes](https://numerical.recipes/)の`locate.f90`をMatlabに翻訳したコードを各自で追加してください。
* `8_3_Khan_Thomas`：8.3節の結果を再現するためのコード
* `8_4_HANK`：テキスト内の結果を再現するためのコード

#### Fortran
* `8_2_Krusell_Smith`：クルセル・スミス・モデル再現用コード
  * 内挿法(Interpolation)に関しては[Numerical Recipes](https://numerical.recipes/)のサブルーチンを利用しています。そのため、ダウンロードしただけではコンパイルは通りません。下記のサブルーチンについて、各自で補完してください：`locate.f90`, `polint.f90`, `splint.f90`, `spline.f90`, `tridag_ser.f90`, `assert_eq2.f90`, `assert_eq3.f90`, `assert_eq4.f90`, `assert_eqn.f90`, `iminloc.f90`, `nrerror.f90`
