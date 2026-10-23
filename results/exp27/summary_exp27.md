# Exp27 — Vakum-Metrik Roket: enerjiyi koruyan roket denklemi

Tarih: 2026-10-16. Kod: `experiments/exp27_vacuum_metric_rocket/vmr27.py`. Komutlar:
- `derive`, `derive_exact`: roket denklemi;
- `casimir`: statik Casimir;
- `dce`: dinamik Casimir;
- `metric`: metrik ayağı;
- `energy`: enerji ayağı;
- `fig`: şekil.

Öngörüler hesaptan önce yazıldı: `predictions_exp27.json`, 14:04. Şekil: `fig_exp27.png`.

**Görev (kullanıcı, /loop).** Eldeki veriyle bir roket denklemi kurmak. Enerjisi bölünen ve dönüştürülebilen hücrelerden gelecek;
itkisi programlanabilir evrensel atomların başka atomları yönetmesinden ve kuantum vakum kontrolünden doğacak. Kısa adıyla
Vakum-Metrik Roket. Üç ayak:
1. enerji geri dönüşümü, η → 1;
2. vakum kontrolü, Casimir ve dinamik Casimir itkisi;
3. programlanabilir atomlarla metrik (uzay-zaman) mühendisliği.

## Kısa sonuç

- **Denklem (türetildi, sayısal olarak sınandı):** β = tanh(η·ln(M₀/M₁)), η = η_dön · η_yön.
- **Kapalı bir araç için:** itki/güç ≤ 1/c = 3,34 nN/W. Vakum ve metrik mekanizmaları bu sınırı aşamıyor.
- **Vakum ayağı:**
  - statik Casimir çevrimi tam olarak sıfır iş yapıyor;
  - dinamik Casimir itkisi en iyi ölçülmüş kaynakta bile ≤ 2·10⁻²¹ N.
- **Metrik ayağı:** warp için gereken negatif enerji, elde edilebilir Casimir enerjisinden 39–71 mertebe büyük. Ayrıca her warp
  sürücüsü yine de bir roket denklemine uyuyor.
- **Doğru kapı:**
  1. η → 1 foton roketi: hücre geri dönüşümü kayıpları azaltıyor, programlanabilir atomlar ışınımı yönlendiriyor (analog metrik);
  2. dışarıdan ışınlanan güç (ışık yelkeni). Yerleşik enerji sınırını aşan tek yol bu.

## 1. Denklem

Kapalı sistem (araç ve yaydığı alanlar) için enerji ve momentum korunur. Aracın anlık durgun çerçevesinde:
- harcanan enerji dE = −c² dM;
- bunun η kesri tek yönde ışınım olarak çıkıyor, momentumu η dE/c;
- kalanı yönsüz çıkıyor, momentumu 0.

Buradan şunlar çıkıyor:

  M dv′ = η dE / c   →   **β = tanh(η ln(M₀/M₁))**;   göreli olmayan sınırda Δv = η c ln(M₀/M₁)

- **İtki sınırı.** Kapalı araç için F ≤ η P/c ≤ P/c (3,34 nN/W). Dışarıdan ışın alan yelken için F = (1+R) P_ışın/c.
- **Enerji kaynağı.** M₀/M₁ = 1/(1 − f ε/c²); f kaynağın kütle payı, ε özgül enerjisi.
- **η = 1** Ackeret foton roketini verir.

**Sayısal sınama (`derive_exact.json`).** Her adımda durgun çerçevede 4-momentumu tam koruyan bir integratör kullanıldı (kalan aracın
enerjisi M − dE, momentumu ηdE), ardından laboratuvar çerçevesine Lorentz dönüşümü yapıldı.
- Toplam enerji ve momentum 10⁻¹² içinde korunuyor.
- η = 1'de kapalı biçimle fark 10⁻¹³.
- η < 1'de hata adım sayısıyla 1/n hızında küçülüyor (20 bin adımda 3·10⁻⁵, 200 bin adımda 3·10⁻⁶). Formül süreklilik limitinde
  doğru.
- İlk kayıtlı koşu (`derive.json`, birinci mertebe Euler) korunumu 10⁻⁵ bozdu. Bunun sebebi yöntemdi, fizik değil.

## 2. Üç ayak

### Enerji geri dönüşümü (η → 1)

Kullanılan hücre ve atom verisi:
- Exp23: hücreler arası aktarım verimi 0,985–1,023; mimari kayıp (yük direnci) %75–83;
- Exp25: atom başına yük çevirme 2,8 meV;
- Exp26: programlanan ²⁹Si başına 63 µJ–1,5 mJ.

Geri dönüşüm kayıpları azaltarak η_dön'ü yükseltir, ama enerji yaratmaz. Kayıp ısı yönlü bir yayıcıdan ışınlanırsa o da itkiye
katılır (ısıl foton itkisi). Kayıp payı 0,8 ve tam yönlenmeyle 1 W girişten 2,7 nN elde edilir.

Δv'yi asıl belirleyen kaynağın kütle-enerji payı. η = 1 ve f = 0,9 için (`energy.json`):

| kaynak | Δv |
|---|---|
| Li-iyon (0,9 MJ/kg) | 2,7 mm/s |
| fisyon (tam yanma) | 246 km/s |
| D-T füzyonu | 1022 km/s |
| antimadde (ideal; ürünlerin yönlendirilmesi gerçekte zor) | 0,98c |

### Vakum kontrolü (Casimir)

- **Statik Casimir.** E = −A π²ħc/(720 a³) bir durum fonksiyonu; bu yüzden (aralık, alan) düzleminde her kapalı çevrimde net iş sıfır.
  - Analitik bacak toplamı 10⁻¹⁶; ilk yamuk kuralı koşusu 4·10⁻⁶ verdi (kuadratür hatası).
  - Levha kuvvetleri iç kuvvetler; aracın net itkisi sıfır.
  - Basınç gerçek: 10 nm'de 1,3·10⁵ Pa.
- **Dinamik Casimir.** Gerçek foton üretir, ama enerjiyi sürücüden alır. Vakum Lorentz değişmez olduğu için momentum vermez;
  F ≤ P_DCE/c.

| kaynak | itki |
|---|---|
| mükemmel ayna (EM): f = ħA q⁽⁵⁾/(30π²c⁴) → P = ħA δq² Ω⁶/(60π²c⁴); 1 GHz, 1 nm | 4,5·10⁻³⁹ N/m² |
| süperiletken devre, Wilson 2011'in ölçtüğü düzey (birim bantta birkaç K) | 7·10⁻²⁴ – 2,3·10⁻²¹ N |
| bölünen hücreler anahtarlanan plazma aynası olarak, ms anahtarlama | 7·10⁻⁷⁸ N |
| karşılaştırma: 1 W foton roketi | 3,3·10⁻⁹ N |

- **Ölçüm.** EmDrive türü "vakum itkisi" iddiaları Tajmar 2021'de sıfır çıktı: olası itki, aynı güçteki klasik ışınım kuvvetinin
  altında.

### Metrik mühendisliği (programlanabilir atomlar)

- **Momentum korunumu.** Kapalı bir sistemin asimptotik momentumu yalnız yayılan momentumla değişir. Bobrick ve Martire (2021):
  her warp sürücüsü itki gerektirir ve bir roket denklemine uyar.
- **Kütleçekim dalgası itkisi.**
  - 1 kg kütle, 1 m, 1 kHz dönüş: 9·10⁻³⁹ N.
  - Hücredeki elektronlar: 7·10⁻¹⁰⁵ N.
- **Negatif enerji.**
  - 100 m'lik Alcubierre balonu ~6,2·10⁶² kg ister (Pfenning ve Ford 1997 sınırıyla; Van Den Broeck 1999).
  - Van Den Broeck'in geometrisi ~3 güneş kütlesi ister.
  - 1 m³'ü 1 nm'lik Casimir aralıklarıyla doldurmak −4,8·10⁻⁹ kg verir; yani 39–71 mertebe eksik.
  - Programlanabilir atomların enerji ölçeği meV.
- **Atomların gerçek katkısı: "analog metrik".** Işınımı yönlendiren programlanabilir bir yüzey (faz dizisi, metayüzey) η_yön'ü 1'e
  yaklaştırır. Korunum yasalarıyla uyumlu metrik katkısı yalnız bu.

## 3. Öngörüler ve sonuç

| no | sonuç |
|---|---|
| VR1 | kısmen. Korunum ✓ (tam adımlı integratör, 10⁻¹²). Formül süreklilik limitinde ✓ (hata ∝ 1/n). Kayıtlı eşik (< 10⁻⁶) 200 bin adımda η < 1 için tutmadı. İlk Euler koşusu korunumu 10⁻⁵ bozdu |
| VR2 | ✓: analitik toplam 10⁻¹⁶ (ilk yamuk koşusu 4·10⁻⁶, kuadratür hatası) |
| VR3 | ✓: 2,3·10⁻²¹ N; 4,5·10⁻³⁹ N/m²; 7·10⁻⁷⁸ N |
| VR4 | ✓: 39–71 mertebe eksik; kütleçekim dalgası itkisi 9·10⁻³⁹ N |
| VR5 | ✓: η → 1 foton roketi; sınır 1/c; ışın kapısı açık |

## 4. Vakum-Metrik Roket Denklemi (son biçim)

Aracın çerçevesinde:

  M dv′/dτ = [ η_dön η_yön P_iç + (1+R) P_ışın ] / c,      dM/dτ = −P_iç / c²

- **η_dön (enerji ayağı).** İç gücün ışınıma dönüşen payı: hücre geri dönüşümü ve yönlü ısı yayıcısı.
- **η_yön (metrik ayağı).** Programlanabilir yüzeyle yönlenme; "analog metrik".
- **P_DCE (vakum ayağı).** P_iç'in bir parçası; ek momentum vermez.
- **Gerçek metrik terimi.** Kapalı sistemde sıfır.

Dış ışın yoksa: β = tanh(η_dön η_yön ln(M₀/M₁)). Sınırlar: F/P_iç ≤ 1/c ve F/P_ışın ≤ 2/c.

## 5. Doğru kapı

- **Kapalı kapılar:**
  - vakumdan bedava momentum ya da enerji (korunum ve Lorentz değişmezliği; ölçüm: Tajmar 2021);
  - statik Casimir çevrimi (durum fonksiyonu);
  - warp metriği (39–71 mertebe; kuantum eşitsizlikleri; Bobrick ve Martire).
- **Açık kapı 1: η → 1 foton roketi.** Hücre geri dönüşümü (η_dön), programlanabilir yönlendirme (η_yön) ve yüksek özgül enerjili
  kaynak birleşiyor. Sınır Δv = c·tanh(ln(M₀/M₁)).
- **Açık kapı 2: dışarıdan ışınlanan güç (ışık yelkeni).** Yerleşik enerji sınırını aşan tek yol. Momentum ve enerji dışarıdan
  geldiği için Δv yerleşik kaynağın kütlesiyle sınırlı değil. F = (1+R)P/c; 100 GW için 667 N.
- **Bizim parçalarımızın yeri.**
  - Hücreler: güç işleme ve kayıp ısının yönlü ışınlanması.
  - Programlanabilir atomlar: ışınımın yönlendirilmesi.

  Hiçbiri itki/güç sınırını (1/c) aşmıyor.

## Sınırlar

- Kuram ve literatür. Vakum ve metrik ayakları için ölçülmüş itki yok (EmDrive sıfır).
- Dinamik Casimir formülü düzlem mükemmel ayna için; gerçek malzeme daha zayıf.
- Antimadde değeri ideal.
- Warp sayıları literatürden; farklı geometriler farklı sayılar verir, ama hepsi elde edilebilir negatif enerjinin çok üstünde.

## Kaynaklar (bu oturumda okundu)

- Dalvit, Maia Neto, Mazzitelli, "Fluctuations, dissipation and the dynamical Casimir effect", arXiv 1006.4790: 1B'de ħq‴/(6πc²);
  3B skaler alanda ħA q⁽⁵⁾/(360π²c⁴); EM alanda TM katkısı TE'nin 11 katı.
- Wilson vd., Nature 479, 376 (2011): sürüş 11,30 GHz; üretilen fotonlar gürültüyü ~2 katına çıkarıyor (birkaç K).
- Tajmar vd., CEAS Space J. 14, 31 (2021): EmDrive itkisi yok.
- Bobrick ve Martire, CQG 38, 105009 (2021): warp sürücüleri itki ister.
- Pfenning ve Ford, CQG 14 (1997): kuantum eşitsizliği; Van Den Broeck 1999: ~6,2·10⁶² kg, değiştirilmiş geometride birkaç güneş
  kütlesi.

## Dosyalar

`predictions_exp27.json`, `derive.json`, `derive_exact.json`, `casimir.json`, `dce.json`, `metric.json`, `energy.json`,
`fig_exp27.png`.
