## 北尾早霧・砂川武貴・山田知明『定量的マクロ経済学と数値計算』日本評論社

### 第8章：定量的マクロ経済学のフロンティアを覗いてみる

本リポジトリは、北尾早霧・砂川武貴・山田知明『[定量的マクロ経済学と数値計算](https://www.nippyo.co.jp/shop/book/9287.html)』（日本評論社）の**第8章**に対応するサポートコードを収録しています。

---

### フォルダ構成

```
chapter8/
├── Julia/
│   ├── 8_2_Krusell_Smith/
│   │   ├── KS_main.ipynb        # メイン（全体の計算を統括）
│   │   ├── KS_household.jl      # 家計の最適化問題（EGM による政策関数）
│   │   ├── KS_markov.jl         # マルコフ連鎖の生成（雇用・TFP の遷移確率）
│   │   ├── KS_simulation.jl     # シミュレーション（集計資本の経路を生成）
│   │   └── KS_OLS.jl            # 集計法則の OLS 推定
│   └── 8_4_HANK/
│       ├── Solve_SS_jl.ipynb    # 8.4節：HANK モデルの定常状態の計算
│       └── Solve_dyn_jl.ipynb   # 8.4節：HANK モデルの動学（IRF の計算）
├── MATLAB/
│   ├── 8_2_Krusell_Smith/
│   │   ├── main_KS.m            # メイン
│   │   ├── end_grid_method.m    # EGM による家計の政策関数の計算
│   │   ├── RHS_Euler.m          # オイラー方程式の右辺の計算
│   │   ├── calibration.m        # カリブレーション
│   │   ├── gen_prob_matrix.m    # 遷移確率行列の生成
│   │   ├── law_of_motion_sim.m  # 集計資本の運動方程式（シミュレーション）
│   │   ├── markov_sim.m         # マルコフ連鎖のシミュレーション
│   │   ├── approx_agg.m         # 集計法則の近似（線形回帰）
│   │   └── regress_KS.m         # 集計法則の OLS 推定
│   ├── 8_3_Khan_Thomas/
│   │   ├── SS/                  # 定常状態の計算（スプライン近似）
│   │   └── KS/                  # Khan-Thomas モデルの完全解（KS アルゴリズム）
│   └── 8_4_HANK/
│       ├── Solve_SS.m           # HANK モデルの定常状態の計算
│       ├── HH_opt_EGM.m         # 家計の最適化（EGM）
│       ├── HH_dist.m            # 家計の定常分布
│       └── Rouwenhorst.m        # Rouwenhorst 法による所得離散化
└── Fortran/
    └── 8_2_Krusell_Smith/
        ├── main.f90             # メイン
        ├── mod_EGM.f90          # EGM モジュール
        ├── mod_calibration.f90  # カリブレーション
        ├── mod_law_of_motion.f90 # 集計資本の運動方程式
        ├── mod_trans_probability.f90 # 遷移確率
        └── （その他モジュール）
```

---

### 各言語のコードについて

#### Julia

* 各節の番号に対応したフォルダにコードが格納されています。
* `8_2_Krusell_Smith`：RAの小野泰輝さんに作成していただいたJulia用再現コードです。もちろん残りうるあらゆる間違いは、すべて著者たちの責任です。
* `8_4_HANK`：8.4節の結果を再現するためのコードです。定常状態の計算（`Solve_SS_jl.ipynb`）と動学的分析・IRF の計算（`Solve_dyn_jl.ipynb`）に分かれています。

#### MATLAB

* 各節の番号に対応したフォルダにコードが格納されています。
* `8_2_Krusell_Smith`：`locate.m` が含まれていないため、途中でエラーが出ます。[Numerical Recipes](https://numerical.recipes/) の `locate.f90` を MATLAB に翻訳したコードを各自で追加してください。
* `8_3_Khan_Thomas`：8.3節（Khan–Thomas モデル）の結果を再現するためのコードです。定常状態の計算（`SS/`）とクルセル・スミス・アルゴリズムによる完全解（`KS/`）の2フォルダに分かれています。
* `8_4_HANK`：HANK（Heterogeneous Agent New Keynesian）モデルの定常状態を再現するためのコードです。

#### Fortran

* `8_2_Krusell_Smith`：クルセル・スミス・モデルの再現用コードです。
  * 内挿法（Interpolation）に関しては [Numerical Recipes](https://numerical.recipes/) のサブルーチンを利用しています。そのため、ダウンロードしただけではコンパイルは通りません。下記のサブルーチンについて、各自で補完してください：`locate.f90`, `polint.f90`, `splint.f90`, `spline.f90`, `tridag_ser.f90`, `assert_eq2.f90`, `assert_eq3.f90`, `assert_eq4.f90`, `assert_eqn.f90`, `iminloc.f90`, `nrerror.f90`

---

### コードの内容

#### `KS_main.ipynb`（Julia）/ `main_KS.m`（MATLAB）　―　8.2節：Krusell–Smith モデル

集計的な不確実性（TFP ショック）と個人的な不確実性（雇用ショック）が共存する Krusell–Smith (1998) モデルを実装します。家計の資産分布を状態変数として扱うことで、一般均衡を解きます。

- カリブレーション：CRRA パラメータ、割引因子、資本減耗率、TFP と雇用のマルコフ連鎖
- 家計の政策関数の計算（EGM を使用）：所与の集計資本の運動方程式のもとで個人の資産選択を最適化
- 集計シミュレーション：マルコフ連鎖に沿って TFP・雇用を発生させ、集計資本の時系列を生成
- 集計法則の OLS 推定：$\log K' = a_0 + a_1 \log K$ の形で集計資本の対数線形則を推定
- Den Haan 精度テスト：モデルの精度を評価する指標の計算

#### `Solve_SS_jl.ipynb`　―　8.4節：HANK 定常状態

Heterogeneous Agent New Keynesian（HANK）モデルの定常均衡を計算します。標準的な Bewley–Aiyagari 型の家計問題と NK のフィリップス曲線・IS 曲線を組み合わせます。

- EGM による家計の最適貯蓄政策関数の計算
- 資産と所得の定常分布の計算
- 財市場・資産市場の清算条件を満たす定常均衡の探索

#### `Solve_dyn_jl.ipynb`　―　8.4節：HANK 動学・インパルス応答

金融政策ショックなどに対するインパルス応答関数（IRF）を計算します。ヤコビアン行列を使った線形化手法（Sequence Space Jacobian など）を実装しています。

- 定常状態周りの線形化
- 金融政策ショック（テイラー則の外れ）に対するインパルス応答の計算
- 消費・資産・利子率・インフレ率の動学プロファイルのプロット
- RANK（代表的家計 NK モデル）との比較

---

### 注意

* Juliaではいくつかのパッケージを利用しています（`Interpolations`、`LinearAlgebra`、`Statistics` 等）。もし実行できない場合はREPLで`]`を押したのち、`add パッケージ名`を実行してインストールしてください。
* 図をプロットする際に日本語が文字化けする可能性があります。数値計算そのものには無関係ですが、日本語で図を表示したい場合はフォントを追加インストールする必要があります。

---

### セットアップ

Julia と Jupyter Notebook のインストール・環境設定については、[インストールと環境構築のガイド](https://quant-macro-book.github.io/setup/) を参照してください。
