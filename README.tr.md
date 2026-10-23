# Güç ve Hücre

**Standart CMOS'ta programlanabilir, kendiliğinden bölünen akım hücreleri: nasıl bölündükleri, nasıl programlandıkları, birbirleriyle nasıl güç alışverişi yaptıkları, fikrin tek tek atomlara kadar ne kadar inebildiği ve denetlenebilir enerji dönüşümüne nasıl büyüyebileceği. Önceden kaydedilmiş tahminli, yayımlanmış ölçümlerle kapatılmış bir model ve benzetim çalışması.**

Araştırma: **counfield** — <https://counfield.com/> · English: [README.md](README.md)

> **Durum (Ekim 2026).**
> - Bu depodaki her şey **model ve benzetimdir**. Proje hiçbir şey üretmedi ya da ölçmedi. Ölçülmüş veri yalnız üçüncü taraftır: modellenen sürecin
>   herkese açık kriyojenik karakterizasyon verisi, yayımlanmış bir silikon ESD akım filamanı veri seti (Mamanee 2012) ve literatürdeki sonuçlar.
> - Her deneyin, koşudan **önce** yazılmış bir `predictions_*.json` dosyası vardır. Tutmayan tahminler saklanır ve raporlanır. Bazı raporlar bir
>   düzeltme kutusuyla başlar; kutu özgün metinle çelişirse kutu geçerlidir.
> - Deney başına raporlar Türkçedir (`results/expN/summary_*.md`). Bu README Türkçe özettir; İngilizce özet `README.md`'dedir.
> - Exp27–30, 34 ve 39–41 (Bölüm 9'daki gelecek hattı) **henüz bu depoda değildir**; birkaç rapor onlara atıf yapar.

---

## Neden "güç ve hücre"

Burada **hücre**, akım taşıyan kendini düzenleyen bir kanaldır. Taşıdığı akımla büyür, çok genişleyince bölünür ve kızları yine bölünür.
**Güç**, hücrelerin taşıdığı ve paylaştığı şeydir: bir hücre öldüğünde komşuları akımını devralır; sürüş değiştiğinde hücreler öngörülebilir bir
biçimde doğar, bölünür ya da yok olur.

Çalışma tek bir sorudan başladı: sıradan bir çip üzerindeki yapı, sabit bir elektriksel sürüş altında kendiliğinden büyüyüp bölünebilir mi? Dört
soruya dönüştü:

1. Hücrenin bölündüğü boyutu ne belirliyor ve yanıt, uydurma sabitsiz bir yasa olarak yazılabilir mi?
2. Doğum, sayı, hayatta kalma ve bölünme elektriksel sürüşle **programlanabilir** mi?
3. Hücre, gücü taşıyan ve geri dönüştüren bir birim olarak ne değer taşıyor ve programlanabilirlik tek tek atomlara kadar ne kadar itilebilir?
4. Denetlenebilir enerji dönüşümünün hangi engellerini programlanabilirlik kaldırabilir? (Bölüm 9, yayımlanacak.)

## Öne çıkanlar

| Sonuç | Sayı | Deneyler | Kanıt |
|---|---|---|---|
| Uydurma sabitsiz kapalı biçimli bölünme yasası, koşulardan önce kaydedilmiş tahminlerle sınandı | 22 bölünmede öngörülen / gözlenen drenaj düşüşü payı **1,001 ± 0,069**; bağımsız yük çözücüsüyle yeniden koşulunca 1,029 ± 0,071 (n = 20) | 17, 20 | model, önceden kayıtlı |
| 10 µm üzerinde +30 mV'luk tohumla seçilen doğum yeri | **10 koşunun 10'u**; doğum yeri kuralı **80** test koşusunun **77**'sinde doğru | 18 | model |
| Ani sürüş basamağının yüksekliğiyle belirlenen hücre sayısı | 1, 2, 3–4 ve 4 hücre; ani zorla küçültme bütün hücreleri aynı anda böler (**3 → 7**) | 19, 20 | model |
| Bir hücre öldürüldüğünde komşulara geçen güç | **%98,5–102,3** (geri dönüşüm verimi 0,985–1,023) | 23 | model |
| Herkese açık 4 K karakterizasyon verisinde yüzer gövde kink'i | n-kanallı **10 aygıtın 10'unda** var (bir aygıtta oran 2,8; 300 K'de 1,15) | 32, 33 | **ölçülmüş (üçüncü taraf)** |
| Tek filamanın akım sınırı | öngörülen **0,889 A**, ölçülen 0,9 A | 31 | **ölçülmüş (üçüncü taraf)** |
| Çok atomlu desen yazma | yalnız küresel sürüş: 7 ALT deseninden 0'ı ve 7 RAND deseninden 0'ı; kapalı çevrim yaz–oku–düzelt: **21/21** | 35, 36 | model |
| Bağımsız ölçümlerle kapanış | 28 sonuç karşılaştırıldı: **16 uyumlu, 6 karışık, 3 uyumsuz**, 3 yöntem | 29–32 | karışık |
| Parametresiz itki yasası, yayımlanmış bir iyon itkici test kümesiyle | **40 işletme noktasının 40'ı** 2σ içinde, RMS **%1,45** | 29 (gelecek hattı) | **yayımlanmış veri** |
| Oda sıcaklığında bölünme | **21 koşuda 0 bölünme** | 21 | model, **olumsuz** |

---

## 1. Soru ve model

### 1.1 Sorunun arkasındaki iki gözlem

- **(i) Kararlı, boyutlu filamanlar.** Bistabil (S biçimli) akım–gerilim karakteristiklerinin yarı iletkenlerde akım filamanları oluşturduğu bilinir
  (Schöll; Hüpper vd. 1993). Yük direnci gibi *küresel* bir kısıt, genişliği akımla büyüyen tek bir filamanı seçer.
- **(ii) Çoğalan lekeler.** *Yerel* bir inhibitörlü tepkime–difüzyon sistemleri kendini çoğaltan lekeler verir (Pearson 1993).

Sorumuz: standart bir CMOS yerleşimi yerel inhibitörü sağlayabilir mi, öyle ki silikondaki akım filamanı bölünen bir hücre gibi davransın? Öyleyse
bölünme boyutunu ne belirliyor?

### 1.2 Aygıt modeli

- **Aygıt.** Standart açık 130 nm CMOS sürecinin geniş bir n-kanallı transistörü; döküm tarafının kompakt model kartından (BSIM4; ngspice ile drenaj
  akımı, çarpmalı iyonlaşma akımı ve eşik gerilimi tablolarına benzetilmiş) kurulan ve açık, **taranmış** kriyojenik varsayımlarla genişletilen
  indirgenmiş bir model: eşik altı eğim, çarpmalı iyonlaşma çarpanı, gövde diyotu başlangıcı, donmuş kuyu sızıntısı ve yanal gövde tabaka direnci.
- **Yerel bistabilite.** Genişliği $dy$ olan her bölümün gövde potansiyeli $V_b$,
  $C_b\,\dot V_b = \kappa\,\partial_y^2 V_b + F_{ii}\,M(V_D)\,j_d(V_b) - j_0 e^{(V_b - V_\text{on})/V_s} - g_\text{sub}V_b$
  denklemini sağlar. Yerel drenaj gerilimi $V_D$'nin bir aralığında bunun bir düşük ve bir mandallanmış durumu vardır: kanal genişliğinin yalnız bir
  kısmında ileten yüzer gövdeli bir mandal, yani **akım filamanı**.
- **İnhibitör.** Drenaj dirençli bir raydır ($r_D$, Ω/µm). Rayın her noktası, dağıtık bir balastla ($r_v$, Ω·µm) ortak bir düğüme $V_c$ bağlıdır.
  $u = V_c - V_D$ düşüşü $u'' = u/\lambda^2 - r_D\,j_d(y)$, $\lambda = \sqrt{r_v/r_D}$ denklemini sağlar.
- **Başvuru kümesi.** SS = 20 mV/dek, $F_{ii}$ = 10, kapı aşırı sürüşü −0,1 V, $V_\text{on}$ = 1,05 V, $V_s$ = 5 mV,
  $g_\text{sub}$ = $10^{-13}$ A/(V·µm); evrenselliği iki küme daha (U1, U2) sınar.

---

## 2. Bölünen filaman (Exp12–13, 16–17, 20)

- Eşiğin hemen altında kutuplanmış geniş bir transistör akım filamanı oluşturur. Drenaj bölümlenir ve her bölüm kendi balast direnciyle besleme
  gerilimine bağlanırsa filaman **akımla büyür, ikiye bölünür ve kızları yine bölünür**.
- **Bölünme yasası.** Filaman, merkezindeki drenaj gerilimi iletken durumun var olabildiği en düşük gerilim olan $V_\text{down}$'a düştüğünde
  bölünür. Yalıtılmış bir filaman için kapalı biçim $w_s = 2\lambda\ln[\lambda/(\lambda-\lambda_c)]$ elde edilir. Sabitleri aygıt modelinden gelir;
  hiçbir şey uydurulmaz.
- **Önceden kayıtlı sınama.** Temiz protokolde (yük $R_L$ = 5 kΩ) gözlenen 22 bölünmede merkezden kenara drenaj düşüşü öngörülen paya
  **1,001 ± 0,069** çarpanı içinde eşit (aralık 0,89–1,11). Çok farklı iki aygıt parametresi kümesinde 0,98 ± 0,13 (9 bölünme). Bütün protokollerde 63
  bölünme çözümlendi. 14 temiz koşunun ilk bölünme genişlikleri kapalı biçimi %9 RMS hatayla izler.
- **Bağımsız yeniden sınama (Exp20).** Özgün çözücünün açık yük doğrusunun, ani bir akım fazlasından sonra dengesiz bir durumda kalabildiği bulundu.
  Örtük yük doğrusuyla Exp17'nin 40 koşusu yinelendi:

  | büyüklük | özgün çözücü | örtük çözücü |
  |---|---|---|
  | mekanizma oranı, temiz | 1,001 ± 0,069 (n = 22) | 1,029 ± 0,071 (n = 20) |
  | mekanizma oranı, sert yük | 1,038 ± 0,062 (n = 22) | 1,003 ± 0,057 (n = 24) |
  | mekanizma oranı, evrensellik kümeleri | 0,979 ± 0,133 (n = 9) | 0,973 ± 0,116 (n = 9) |
  | ilk bölünme genişlikleri, temiz | n = 14, ort. +%1,6, RMS %9,1 | n = 13, ort. +%3,2, RMS %8,6 |

  Yayımlanmış 14 ilk bölünmenin 13'ü yeniden üretiliyor. 5,38 V'ta bölünen (%+20) bir koşu artık 5,6 V'a kadar bölünmüyor ve çok düşük sızıntılı bir
  koşu için "bölünme yok" ifadesi geri çekildi (175 µm'de bölünüyor, öngörünün yaklaşık iki katı).
- **Bir tahmin tutmadı.** λ = 36 µm'de, yani $\lambda_c$'nin %1 üstünde, filaman bölünmedi.
- **Denklem keşfi (Exp16).** Seyrek regresyon hattı önce bilinen denklemlerde doğrulandı: temiz veriden bir tepkime–difüzyon sistemini tam olarak
  buluyor, güçlü formda %1 gürültüde çöküyor ve çarpmalı iyonlaşma yasasını 90 eğrinin 90'ında buluyor. Filaman verisinde $w = I/j_\text{on}$
  genişlik yasasını seçiyor ($j$ ≈ 2,55 ± 0,03 µA/µm, R² ≥ 0,9988). Tek serbest sabitli bir uydurma kuramsal değerin 1,01 katını veriyor. Kör sembolik
  regresyon benzer sıkı bir yasa bulamadı. İki sınama başarısız oldu ve raporlandı: bir sonlu elemanlar küre kapasitansı yakınsaması ve verinin karara
  bağlayamadığı bir ölçekleme sınaması.

## 3. Boyut bölünmeyi belirler (Exp11, 14, 15)

- **Bir boyut penceresi.** Kimyasal olarak aktif damlaların süreklilik modellerinde bölünme yalnız bir boyut penceresi içinde olur. Pencerenin altında
  (d ≈ 2–2,2) damla şerit olur; üstünde damla içten boşalıp sıvı kabuğa dönüşür. Pencere, iç tüketim hızı $\Delta\nu$ arttıkça yukarı kayar ve
  genişler: $\Delta\nu$ = −0,004'te d ≈ 2,25–2,5; −0,001'de 2,3–3,5; +0,004'te 2,5–4 ve $\Delta\nu \ge +0{,}012$ için 6B damlalar iki kompakt kıza bölünür.
- **Gerçek bir 6B koşusu.** Bir Cahn–Hilliard aktif damlası 6B hiperküpte (16⁶ = 16,8 milyon nokta) koşuldu. O(5) simetrik bir indirgeme bunu t = 40'a
  kadar ≤ %2,5 farkla yeniden üretiyor ve yüksek boyutlu bölünme ve kabuk sonuçları bu doğrulanmış indirgemeden geliyor.
- **3B'ye karşı 2B.** 3B aktif damla taranan bütün sürüş aralığında iki eşit kıza bölünür ve bölünme yarıçapı sürüşle sürekli değişir. 2B'de aynı çözücü
  solucan ve labirent verir (ikili bölünme 30 koşuda 0; kalıcı bölünme süperkapasitör türü modelde 28'de 0'a karşı 30'da 2, Fisher p = 0,49, anlamsız).
- **Otokatalitik lekeler.** Bir Gray–Scott taramasında (F, k) noktalarının **42'sinin 5'inde** çoğalan lekeler sürer; 6'sının 6'sı %5 gürültüye dayanır ve
  yeni lekeler bölünmeyle doğar (arka plandan doğum yok).
- **Transistör kanalı 2B ve kapılıdır.** Orada süperkapasitör türü ve aktif damla mekanizmaları hücre benzeri bölünme vermez; yalnız otokatalitik
  sistemler aday kalır. Silikonda bölünen geometri tek boyutludur: 2B dizilerde akım lekesi bölünmek yerine solucana uzar.
- **Faz uzayı delikleri (Exp15).** Çarpışmasız bir plazmada elektron delikleri 1B'de birleşir, 2B ve 3B'de parçalanıp yok olur. Hiçbir boyutta kalıcı
  bölünme bulunmadı; bu klasik sonuçlarla uyumludur (Morse ve Nielson 1969).

## 4. Hücreleri programlamak (Exp18–21)

| Kol | Sonuç |
|---|---|
| Hücrenin **nerede** doğduğu | İlk hücre aşırı sürüşü en büyük tek hücrede doğar (kalibrasyon koşularının 19/20'si). Kural doğum yerini **77/80** test koşusunda doğru öngördü. 10 µm üzerinde +30 mV'luk bir tohum, ilk hücrenin yalnız orada doğmasını **10/10** koşuda sağladı; gereken pay, hücre başına rastgele eşik yayılımının ($\sigma$ = 12 mV / 0,5 µm) yaklaşık 2,5 katıdır. |
| **Kaç** hücre doğduğu | $V_s$ yüksekliğindeki ani bir sürüş basamağı sayıyı belirler: $V_s \le 1{,}0$ V 1 hücre, 1,05–1,1 V tam 2 (10 alanın 10'unda), 1,2 V 3–4, 1,4 V dört tohumun hepsi (10/10). Tohum dışında hiç hücre doğmadı. Yavaş taramada hep bir tohum kazanır (halka cihazda en güçlü giriş 19/20). |
| **Hangilerinin** yaşadığı | İki düğme: darbe yüksekliği doğumları, tutma gerilimi yaşayanları (bir taşıma kapasitesi) belirler. Çip üstü ortak düğüm kapasiteleri ($R_L C \le 5$ ms) sonucu değiştirmez. |
| **Hasar ve yenilenme** | Yavaş ölen bir hücre akımını komşulara verir, komşular genişler (3,3 V'ta 3 → 2 hücre); yeni hücre doğmaz (16 koşunun 0'ında). Ani ölüm hasar kenarlarında 2–3 yeni hücre doğurabilir (3 → 4–5). Komşunun merkezini 150 mV zayıflatmak onu hemen böler (3 → 2 → 3). 5,0 V'ta genişleyen komşu kaybolan hücreyi kendiliğinden yeniler; bölünme düzeni oranı 0,93. |
| **Hücreyi taşımak** | Hareketli bir kapı deseni (önde +40 mV, arkada −40 mV) hücreyi aygıt ucuna 5 adımda (3,3 V) ya da 9 adımda (5,0 V) taşır. |
| **Aynı anda bölmek** | Ani zorla küçültme (bir kapı çukurunu bırakmak ya da en az 0,2 V'luk bir $V_{DD}$ düşüşü) hücreleri eşzamanlı böler: 5,0 V, −60 mV'ta **3 → 7**. |

**Yol boyunca bulunan sınırlar.**
- Sürüş, doğumda *k kazananlı* değil *tek kazananlı* bir seçicidir: iki eşit tohum 10 koşunun yalnız 2'sinde birlikte doğdu ve bir basamağın
  kazananları en güçlü k tohum değildir. Sonradan bulunan bir öngörücü (eşik altı toplam ağırlık) daha iyidir ama her yeni alanda kayıtlı eşiği tutmadı.
- **300 K'de 21 koşuda bölünme olmadı (Exp21).** Bistabilite zayıf (pencere 43–68 mV; 4 K'de yaklaşık 770 mV; açık/kapalı akım oranı 5–12; 4 K'de yaklaşık
  $10^4$); kuyu donmuyor, bu yüzden 600 µm'lik gövdenin tamamı tek düğüm gibi davranıyor; ve aygıt uyumsuzluğu ince pencereyi yiyor. 4 K'de bölünme geniş
  bir bistabil pencere ve bağımsız mikro hücreler ister.

## 5. Hücreler arasında enerji (Exp22–23)

- **Alt bantlar.** Evirtim tabakası alt bantları sürecin her transistöründe var (300 K'de elektronların %47–72'si en alt bantta). Coulomb blokajı
  ≤ 0,4–0,8 K ister; filamanların kendisi klasiktir. Bu aygıt sınıfındaki dissipatif bir hücre kuantum noktasının gerektirdiği boyuta inemez (4 K'de
  kuantum hedefi 6–20 nm hücre yarıçapı ister; en iyimser dissipatif hücre 2B tabakada yaklaşık 123 nm'dir).
- **Güç alışverişi.** Bir hücre öldürüldüğünde gücünün **%98,5–102,3'ü** yaşayan hücrelerde yeniden ortaya çıkar (3,3 V ve 5,0 V'ta yavaş ve ani
  ölümde geri dönüşüm verimi η = 0,985–1,023; en kötü durum, 3,3 V'ta yavaş ölüm, %1,5 kaybeder). Her komşunun aldığı pay yakınlığa bağlıdır (iki komşuda
  yaklaşık %73 / %27) ve bir kapı ön gerilimi bunu yalnız ±%4 değiştirir.
- **Kenetlenme ve açık bir uyarı.** Hücre gerilimi $V_M$'de kenetlenir, yani hücre enerjiyi yoğunlaştırmaz. Geri dönüşüm neredeyse kayıpsızdır,
  **ama mimari değildir**: kaynak gücünün %74,6'sı (3,3 V) ve %83'ü (5,0 V) yük direncinde harcanır ve hücrelere yalnız %24 ve %16 ulaşır. Bir hücre
  neredeyse hiç enerji de depolamaz (gövde enerjisi $10^{-14}$ J mertebesinde): bir *akış* yapısıdır, enerji deposu değil. Bu yük kaybını azaltmak
  açık bir tasarım sorunudur (Bölüm 9).

## 6. Hücrelerden atomlara (Exp24–26, 35–36)

**6.1 Hücre atom boyuna inemez (olumsuz).** Üç taban türetildi ve sınandı: donmuş gövde bağlanma uzunluğu (bu modelde yaklaşık 2,2 µm), tükenme genişliği
(yaklaşık 83 nm) ve en az 36 nm'lik 10 yıllık kararlılık genişliği. Programlanabilirlik bağlanma uzunluğunun altında kayboldu. "Programlanabilen hücre
atom kadar küçülür" önermesi elektronik filaman için yanlıştır; oksitteki iyonik filaman için doğrudur (kinetik Monte Carlo'da tek atom genişliği), ama
bu farklı bir aygıt sınıfıdır.

**6.2 Hücreye hapsolmuş atomlar 0/1 yapılabilir (model).** Donmuş hücre gövdesi, ya nötr ya iyonize akseptör atomları taşır; her atomun kendi anahtarlama
gerilimi, yani "adresi" vardır. 1 K'de 100 nm'lik bir gövdede her iki yanda en az 8 kT/e paylı yaklaşık 35 atom vardır ve yaklaşık 31'i %88 başarıyla
yazılır. 4 K'de çalışan yöntem, adresleri 4 K'de kalibre edilmiş yaklaşık 60 nm'lik bir gövdedir: yaklaşık 3,3 atom ≈ %99 doğrulukla programlanır.
40 nm'de adresler pencerenin dışına kayar. Kayıtlı ısıl tahminler (ortalama doğruluk ≥ 0,9) tutmadı: 4 K'de 0,82 ve 1 K'de 0,88.

**6.3 Alanlar mantığı programlar (model, 4 K).** Hücre alanı (0,4–1,8 V/µm) bir donör atomunun aşırı ince ve Zeeman enerjilerini kHz–MHz kaydırır; atomu
frekansla seçmeye yeter. Üç doğal akseptör atomundan kurulan bir NAND kapısı, eşikler her giriş durumu için kalibre edilince 28 üçlünün 27'sinde ve 10
adanın 9'unda çalışır. İlk yöntemde bir pencere hatası (28'de 3) vardı ve kayıtta durur.

**6.4 Evrensel atom (model).** Merkezdeki bir donör elektronu yakın çekirdek spinleriyle etkileşir ve onları tam bir kapı setiyle sürebilir. Tek bir RF
bobini elektron + 4 çekirdek GHZ dolanıklığını **0,999** doğrulukla kurar; yaklaşık 2000 çekirdek daha dahil edilince 0,987; en az 10 ms elektron
eşevreliliğiyle yaklaşık 0,91. Yayımlanmış silikon GHZ doğrulukları %92,5–96,2'dir; yani model iyimserdir (Bölüm 7'deki üç uyumsuzluktan biri).

**6.5 Doğal ve saflaştırılmış silisyum.** Doğal silisyumda güçlü bağlı komşuların yalnız %1–3'ü tek tek seçilebilir (spektral kalabalık) ve depolama ölçülmüş
koşullarda çöker. **Her** komşuya erişen evrensel bir programlanabilir atom saflaştırılmış ²⁸Si'de mümkündür (50 ppm'de 6 nm içindeki bütün ²⁹Si'ler
yapılandırmaların %98'inde) ve doğal silisyumda değildir; engel en yakın komşu ²⁹Si–²⁹Si çiftleridir ve saflaştırma onları ortadan kaldırır. Saflaştırılmış
silisyumda ve elmasta yayımlanmış deneyler (Exp26 raporunda anılır) malzeme koşuluyla uyumludur.

**6.6 Kontrol yığını (Exp35–36).**
- **Tek bir küresel sürüş keyfî bir desen yazamaz:** tekdüze desenler 7/7, ALT 0/7, RAND 0/7.
- **Küresel düzeltimli kapalı çevrim yaz–oku–düzelt** hedef desene **21/21** görevde ulaşır (seçilmiş 3–4 atomlu desenler, modal), marjinal 0,992, tam
  desen 0,973. Yalnız açık çevrim adres kalibrasyonu başarısız oldu (zor desenlerde 14'te 2).
- **Bedel arka plan bozulmasıdır:** hedef dışı atomların %6,25'i çevriliyor; kayıtlı sınır %5. Sonuç bütün gövdenin yalıtımı değil, *seçilmiş bir hedef
  altkümesinde* tam kontroldür.
- **Uzamsal çözünürlük belirleyicidir:** fiziksel bir yüzey çekirdeği ve tek atışlık okumayla 6×6 ızgara kayıtlı tam-kontrol kuralını geçer (modal 13/14,
  arka plan %4,70); daha ince adımlı 8×8 ızgara 14/14 verir (arka plan %3,95).

## 7. Ölçülmüş veriyle kapanış (Exp31–33)

- **Silikonda akım filamanları (Mamanee 2012, 300 K).** Sabit filaman sayısında filaman başına akım yoğunluğu sabittir (%0,1 ve %5,2), sayı yasası
  $N \approx I/(0{,}75\,j\,w_s)$ ölçülen aralıkların içine düşer (ölçülen ortalamanın 0,84–0,90'ı) ve tek filaman sınırı ölçülen 0,9 A'e karşı **0,889 A**'dir.
  **Yasa bir noktada tutmuyor:** akım yoğunluğu N = 2'den N = 4'e %22 düşüyor (modelde tek bir $j_\text{on}$ vardı).
- **Modellenen sürecin herkese açık 4 K karakterizasyonu.** Yüzer gövde kink'i 4 K'de **10 aygıtın 10'unda** var (W7/L0,15 aygıtında oran 2,8; 300 K'de
  1,15). Ölçülen çarpmalı iyonlaşma akımı ve transfer eğrisiyle hücre sabitleri $V_M$ = 0,70 V ve $j_\text{on}$ = 64 µA/µm olur (model sabitleri, filamanın
  ölçümü değil). Gövde etkisi Exp12 modelinin varsaydığından 2–3 kat küçüktür (yalnız ters yönde ölçüldü).
- **Ölçülmeyen:** filamanın kendisi ve histerezisi. Herkese açık en geniş aygıt 7 µm'dir, balastsızdır ve tek yönde taranmıştır.
- **Sayım tablosu (Exp32).** Exp0–27'nin 28 sonucunun hepsi bağımsız yayımlanmış bir ölçümle karşılaştırıldı: **16 uyumlu** (ölçümle düzeltilenler dahil),
  **6 karışık ya da kısmen**, **3 uyumsuz** (Exp22'nin alt bant çözücüsü, Exp25'in benzetilmiş GHZ'si ve Exp26'nın benzetilmiş eşevreliliği) ve 3 yöntem deneyi.

## 8. Çalışmayanlar (bilerek saklanıyor)

- Oda sıcaklığında *kuantum* hücre denemesi (Exp0–10) şimdilik başarısız: ölçek, bağlaşım ve tutarlılık gereken rejime ulaşmıyor.
- 300 K'de bölünme yok (Exp21); faz uzayı deliklerinde kalıcı bölünme yok (Exp15); transistör kanalı damlasında 2B'de bölünme yok (Exp11).
- Bölünme yasası filaman sayısı bağımlılığında (Exp31) ve λ = 36 µm öngörüsünde (Exp17) tutmuyor; bir ilk bölünme ve bir "bölünme yok" ifadesi bağımsız
  yeniden sınamadan sonra değişti (Exp20).
- Kaynak gücünün büyük kısmı mevcut mimaride hücrelerin dışında harcanıyor (Exp23).
- Programlanabilir hücre silikonda atom boyuna inemiyor (Exp24); atom programlamanın kayıtlı ısıl doğruluğu tutmadı.
- Üç model sonucu ölçümlerle uyuşmuyor (Exp22, 25, 26); açık çevrim adresleme başarısız oldu ve küresel düzeltimli çevrim kayıtlı arka plan sınırını aşıyor
  (Exp35).
- Denklem keşfi çalışmasının birkaç çözümleme adımı başarısız oldu ya da sınanamadı (Exp16).

## 9. Gelecek: kontrol edilebilir "yakıt hücreleri" (yayımlanacak)

> **Durum.** Exp27–30, 34 ve 39–41, programlanabilir hücreleri **denetlenebilir enerji hücreleri** olarak ele alan ayrı bir hattır ("yakıt hücresi", kontrollü bir
> miktarda güç dönüştürüp salan birimler anlamında; elektrokimyasal hücre değil). Kuram düzeyindedir; yalnız bir kısmı yayımlanmış veriyle kapatıldı. Raporlar
> tahminleri ve verileriyle sonraki bir güncellemede yayımlanacak. Aşağıdaki sayılar bu raporlardan alıntıdır.

**9.1 Fikir.** Hücreler pozitif enerji sağlar; atomlar (iyonlar) momentumu taşır. Hücre, karakteristik akımı 0,2 mA mertebesinde, kendini düzenleyen bir akım
kanalıdır; hücre sayısı sunulan gücü basamaklarla belirler; ölen bir hücre gücünü komşularına verir.

**9.2 Şimdiden kurulmuş olanlar.**
- **Tepkime kütlesiz yollar kapalı.** Kapalı bir araç 1/c = 3,34 nN/W'lık itki/güç oranını aşamaz; vakum dalgalanması yolları en iyi ölçülmüş kaynakta en çok
  2 × 10⁻²¹ N verir; negatif enerjili metrik yolları erişilebilir herhangi bir vakum enerjisinden 39–71 mertebe fazla enerji ister (Exp27, olumsuz sonuç).
- **Atom yolu birim enerjiye çok daha fazla momentum taşır.** Harcanan birim enerji başına atom 2/v_e, foton 1/c momentum taşır; oran 2c/v_e: 60 km/s'de
  $10^4$, 10 km/s'de $6 \times 10^4$ (Exp28).
- **Parametresiz bir itki yasası, $F/P = 2\eta/v_e$ (demet için $F = \gamma J_b\sqrt{2mV/q}$), yayımlanmış veriyle sınandı** (Exp29–30). Yayımlanmış bir iyon
  itkicinin tezgâhta ölçülen **40 işletme noktasının 40'ını** serbest parametre olmadan 2σ içinde öngörür (ortalama sapma +%1,0, RMS **%1,45**; 1000 V üstündeki 33
  noktada %1,18), uçuşta ölçülmüş bir noktayı (92 mN, 2,3 kW) η = 0,61 ile yeniden üretir ve bir sıvı-metal alan emisyonlu itkicinin doğrudan ölçülen itkisini
  %1,6–4,3 içinde verir.
- **Atom düzeyindeki halkalar yayımlanmış deneylerde ölçülmüştür:** istek üzerine tek iyon (%88,5, algılama sınırlı), uçuş süresinin $\sqrt{m}$ ölçeklemesi
  (1,97σ) ve yapay atomda $I = e f$ (0,27 ppm).
- **Model zarfı.** Modelde her hücre 7–26 µN verir (1–5 kV), hücre sayısı itkiyi basamaklarla belirler ve bir hücre ölünce itki en çok %0,8 değişir (Exp28).
  Önceki bütün verilerden enerji bakımından en iyi bir işletme zarfı türetildi (Exp34); tutmayan alt öngörüler dahildir (hesap hatasıyla $10^3$ kat kaçan bir
  kayıtlı aralık raporda düzeltilmiştir).
- **Keşif çalışmaları.** Görev ölçeğinde ölçekleme ve aşırı ölçek inişi çalışmaları (Exp39–41) spekülatiftir ve bu etiketle yayımlanacak.

**9.3 Programlanabilirliğin kaldırabileceği engeller.**

| Denetlenebilir enerji dönüşümünün engeli | Programlanabilir kol | Buradaki kanıt | Durum |
|---|---|---|---|
| Güç ince, geri döndürülebilir basamaklarla kısılamıyor | hücre sayısı basamak yüksekliğiyle, yaşayanlar tutma gerilimiyle belirlenir; zorla küçültme bütün hücreleri böler (3 → 7) | Exp19–20 | model; 300 K'de bölünme yok (Exp21) |
| Tek bir yerel arıza bütün birimi düşürüyor | ölen hücre gücünün %98,5–102,3'ünü komşulara verir; yavaş ölüm yeni hücre doğurmaz (16'da 0); kasıtlı hasar bölünmeyle onarılır (3 → 2 → 3) | Exp18, 23 | model |
| Akım yığılması ve denetimsiz filamanlaşma | küresel kısıt tek filaman seçer; dağıtık balast onu öngörülebilir bir genişlikte böldürür, $w = I/j_\text{on}$ | Exp12–17 | model; ESD filamanları ölçülmüş (Exp31) |
| Etkinliğin başladığı yer rastgele | +30 mV'luk tohum doğum yerini sabitler (10/10; kural 77/80) | Exp18 | model |
| Bir birim kapatılınca enerji kaybolur | hücreler arası geri dönüşüm neredeyse kayıpsız (η 0,985–1,023) ve hücre gerilimi $V_M$'de kenetli; **ama kaynak gücünün %75–83'ü yük direncinde kayboluyor** | Exp23 | model; **açık** |
| Programlanabilirlik geri besleme ister | yaz–oku–düzelt çevrimi hedef desene 21/21 görevde ulaşır; yalnız küresel sürüş ulaşamaz (0/7) | Exp35–36 | model |
| Birim yeterince küçük yapılamıyor | elektronik hücrenin tabanları var (≈ 36 nm kararlılık, ≈ 2,2 µm gövde bağlanması); gövdedeki programlanabilir akseptör atomları daha küçük birimdir | Exp24 | model; tabanlar sınırdır, aşılmadı |
| Malzeme saflığı kontrolü sınırlıyor | bütün komşuların evrensel kontrolü saflaştırılmış malzeme ister | Exp26 | model + anılan ölçümler |
| Oda sıcaklığında çalışma | elektronik filaman 300 K'de bölünmez; oksitteki iyonik filaman aday sınıftır | Exp21, 24 | olumsuz; yalnız benzetim |
| Doğrulama boşluğu | 28 sonuç bağımsız ölçümlerle karşılaştırıldı; itki halkası kapalı (40/40) | Exp29–32 | kısmen kapalı |

**9.4 Hâlâ açık olanlar (ve sıradakiler).**
1. Zincirin ölçülmemiş tek halkası **hücrenin kendisi, kriyojenik sıcaklıkta**: tohumlu, balastlı bir aygıtın akım basamakları.
2. Gerçek bir aygıtta tohumlu doğum yeri sınaması (Exp18 bir modeldir; silikon filamanı ölçümlerinde tohumlu aygıtların doğum yerini sabitlediği bildirilir).
3. Bölüm 5'teki yük kaybını azaltmanın bir yolu, oda sıcaklığında bölünen bir aygıt sınıfı (iyonik filamanlar adaydır) ve Bölüm 6.6 için yeterince hızlı bir kapalı
   çevrim okuma.
4. Bu bölümdeki her sayı yeniden üretilebilsin diye Exp27–30, 34 ve 39–41'in verileriyle birlikte yayımlanması.

## 10. Yöntem: önceden kayıt ve dürüstlük

- **Koşulardan önce tahmin.** Her deney benzetimden önce bir `predictions_*.json` yazar; dosya bir kez yazılır ve sonraki betikler onu ezmemelidir. Tutmayan
  tahminler, tutanların yanında saklanır ve raporlanır.
- **Açıklama.** Bir tahmin kör olmadığında (örneğin kaynaklar ölçütler yazılmadan önce okunduysa) rapor bunu söyler.
- **Düzeltmeler görünür kalır.** Denetimler ve yeniden sınamalar ilgili raporların başına kutu olarak eklenir; kutunun altındaki metin kayıt için korunur.
- **Kapanış.** Her model sonucu, varsa bağımsız bir ölçümle karşılaştırılır; yoksa "ölçüm yok" ya da "yöntem" diye işaretlenir. Karşılaştırmaya uydurma yapılmaz.
- **Olumsuz sonuçlar** bir dipnot değil, sonuç kümesinin parçasıdır (Bölüm 8).

## 11. Depo yapısı ve yeniden üretme

| Yol | İçerik |
|---|---|
| `src/qcell/` | modeller: sonlu elemanlar, Poisson–Schrödinger, faz alanı, aktif madde, aygıt tabloları |
| `experiments/` | deney başına bir klasör (`expN_*`), çalıştırılabilir betikler |
| `results/` | `expN/`: Türkçe özet raporlar, `predictions_*.json` (koşudan önce yazılmış), veri, şekiller |
| `tests/` | 74 pytest testi |
| `data/` | üçüncü taraf aygıt modellerinin ve ölçüm verisinin kayıtları ve lisansları |

```bash
python -m venv .venv && source .venv/bin/activate     # Windows: .venv\Scripts\activate
pip install -r requirements.txt
PYTHONPATH=src pytest -q                              # 74 test
```

Denklem keşfi deneyleri ayrı bir ortamda `requirements-optional.txt` ister. Bazı deneyler `ngspice` (aygıt tabloları; tablolar dahildir) ya da FEniCSx/dolfinx
(isteğe bağlı çapraz kontroller) ister; her deneyin raporu kullandıklarını sayar. Büyük ham diziler saklanmaz; betikler yeniden koşulunca üretilir.

## 12. Sınırlar

- **Model, ölçüm değil.** 4 K aygıt davranışı 300 K kompakt model tablolarından ve taranmış varsayımlardan bir ekstrapolasyondur. Dahil değil: öz-ısınma, tuzak
  kinetiği, parazitik bipolar transistör ve kaynak tarafı direnci.
- **Sürüş ve yük.** Bölünme yasası yarı statik taramayı varsayar. Daha sert bir yük ($R_L$ = 2 kΩ) bölünmeleri yaklaşık %4, düşük uyumsuzlukta %14 geciktirir; çok
  düşük gövde sızıntısı geciktirir ya da bastırır.
- **Saçılma ve eşik.** Tekil bölünmelerde ölçüt ±%7 saçılır (1σ, temiz protokol). Eşiğe yakın yerde öngörülen $w_s$ sabitlere aşırı duyarlıdır.
- **Sınanmamış parçalar.** Periyodik dizi formülü ve filaman sayısı kılavuz çizgisi doğrulanmadı.
- **2B sonuçlar.** 2B koşularda gerilim adımlarının çoğu kararlı durum toleransına ulaşmadı; 2B sonucu (bölünme yok) 1B sonuçlardan zayıf kanıttır. $d \ge 4$ aktif
  damla sonuçları O(d−1) simetrik indirgemeden gelir; yalnız 3B'de ve büyüme sırasında 6B hiperküple doğrulandı.
- **Atom düzeyi sonuçlar modeldir.** İdeal kontrol alanlarını varsayarlar ve kapanış tablosu nerede abarttıklarını gösterir (Bölüm 7).
- **Kuantum değil.** Hücre sonuçları klasik, dissipatif desen oluşumudur.

## 13. Seçilmiş kaynaklar

- A. Hüpper, E. Schöll vd., Phys. Rev. B 47, 15515 (1993); E. Schöll ve D. Drasdo (1990): S biçimli akım–gerilim karakteristikleri ve akım filamanları.
- J. E. Pearson, *Complex patterns in a simple system*, Science 261, 189 (1993).
- D. Zwicker vd., *Growth and division of active droplets provides a model for protocells*, Nat. Phys. 13, 408 (2017).
- W. Mamanee, *Investigation of current filament behavior in smart power ESD devices*, TU Wien tezi (2012).
- B. Morse ve C. Nielson, Phys. Rev. Lett. 23, 1087 (1969): faz uzayı delikleri.
- T. D. Ladd vd., Phys. Rev. B 71, 014401 (2005): silikonda nükleer spin eşevreliliği.
- Exp26 raporunda anılan yayımlanmış atom düzeyi deneyler: Pla vd. (2014), Asaad vd. (2020), Mądzik vd. (2022), Thorvaldson vd. (2025), Edlbauer vd. (2025),
  Bradley vd. (2019).

## 14. Lisans, üçüncü taraf veri ve atıf

MIT lisansı (bkz. `LICENSE`). Aygıt modeli dosyaları ve herkese açık ölçüm kayıtları kendi lisansları altındaki üçüncü taraf malzemedir; bkz. `NOTICE.md`. Yayımlanmış
makalelere atıf yapılır, yeniden dağıtılmaz. Atıf için `CITATION.cff`'e bakın.
