/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hw;

import java.awt.Graphics;
import java.awt.Image;
import java.awt.image.MemoryImageSource;
import java.awt.image.PixelGrabber;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.imageio.ImageIO;
import javax.swing.JFrame;
import javax.swing.JPanel;

/**
 *
 * @author Qbold
 */
public class HW extends JFrame {

    public static void main(String[] args) {
        new HW().setVisible(true);
    }

    private static Image back, fore;
    private final static int ANCHOR_X = 184, ANCHOR_Y = 62;
    private static double eps = 1;

    public HW() {
        setTitle("Poisson");
        setDefaultCloseOperation(EXIT_ON_CLOSE);
        add(new JPanel() {

            @Override
            public void paint(Graphics g) {
                HW.paintPanel(g);
            }

            @Override
            public void update(Graphics g) {
                paint(g);
            }
        });

        loadImage();
    }

    private void loadImage() {
        try {
            back = ImageIO.read(new File("BarackObama.jpg"));
            fore = ImageIO.read(new File("mona.png"));
            fore = getPoisson(fore, back);

            int[] p = getPixels(back);
            for (int i = 0; i < p.length; i++) {
                int r = (p[i] >> 16) & 0xff;
                int g = (p[i] >> 8) & 0xff;
                int b = p[i] & 0xff;

                r /= 1.2;
                g /= 1.2;
                b /= 1.2;

                p[i] = 0xff000000 | (r << 16) | (g << 8) | b;
            }
            back = createImage(new MemoryImageSource(back.getWidth(null), back.getHeight(null), p, 0, back.getWidth(null)));

            setSize(back.getWidth(null), back.getHeight(null));
        } catch (IOException ex) {
            Logger.getLogger(HW.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    private static void paintPanel(Graphics g) {
        if (back != null) {
            g.drawImage(back, 0, 0, null);
            g.drawImage(fore, ANCHOR_X, ANCHOR_Y, null);
        }
    }

    private int[] getPixels(Image w) {
        int[] pix = new int[w.getWidth(null) * w.getHeight(null)];
        PixelGrabber gr = new PixelGrabber(w, 0, 0, w.getWidth(null), w.getHeight(null), pix, 0, w.getWidth(null));
        try {
            gr.grabPixels();
        } catch (InterruptedException ex) {
            Logger.getLogger(HW.class.getName()).log(Level.SEVERE, null, ex);
        }
        return pix;
    }

    private int[] getChannel(int[] pix, int bitshift) {
        int[] bts = new int[pix.length];
        for (int i = 0; i < pix.length; i++) {
            bts[i] = (pix[i] >> bitshift) & 0xff;
        }
        return bts;
    }

    private int[] getH(Image fon, int[] pix, int w, int anch_x, int anch_y) {
        int[] data = new int[pix.length];
        int h = data.length / w;

        int[] p = getPixels(fon); // пиксели фона

        //   p = getGrayScale(p);
        //  mona_lisa = createImage(new MemoryImageSource(fon.getWidth(null), fon.getHeight(null), p, 0, fon.getWidth(null)));
        for (int i = 0; i < pix.length; i++) {
            if (pix[i] > 0) {
                int x = i % w;
                int y = i / w;

                int x_a = x + anch_x;
                int y_a = y + anch_y;

                if (x > 0 && pix[i - 1] == 0 || x + 1 < w && pix[i + 1] == 0 || y > 0 && pix[i - w] == 0 || y + 1 < h && pix[i + w] == 0) {
                    data[i] = p[x_a + y_a * fon.getWidth(null)];
                }
            }
        }

        return data;
    }

    private int[] mergeChannels(int[] a, int[] r, int[] g, int[] b) {
        int[] res = new int[a.length];

        for (int i = 0; i < res.length; i++) {
            res[i] = (a[i] << 24) | (r[i] << 16) | (g[i] << 8) | b[i];
        }

        return res;
    }

    // alpha - маска накладываемой области
    // color - цвет накладываемой области
    // h_alpha - маска границы
    // h_color - цвет границы
    // w - ширина изображения
    private int[] runPoisson(int[] alpha, int[] color, int[] h_alpha, int[] h_color, int w) {

        int h = alpha.length / w;
        int[] res = new int[alpha.length];

        ArrayList<Pixel> pixel = new ArrayList<>(alpha.length);
        int[] index_to_number = new int[alpha.length]; // номер переменной по индексу пиксела в изображении
        int[] number_to_index = new int[alpha.length]; // индексы пикселов в изображении по номеру переменной в системе уравнений
        for (int i = 0; i < alpha.length; i++) {
            if (alpha[i] > 0) {
                index_to_number[i] = pixel.size();
                number_to_index[pixel.size()] = i;
                pixel.add(new Pixel(i % w, i / w));
            } else {
                index_to_number[i] = -1; // если данный пиксел не принадлежит рассматриваемой области
            }
        }
        for (int i = alpha.length - 1; number_to_index[i] == 0; i--) {
            number_to_index[i] = -1;
        }

        int[] grad = getGradient(color, w); // градиент значений текущего канала изображения

        int[][] D = {{1, 0}, {-1, 0}, {0, 1}, {0, -1}}; // массив векторов (dx; dy)

        MatrixLine[] lines = new MatrixLine[pixel.size()]; // строки матрицы коэффициентов системы линейных уравнений
        double b[] = new double[lines.length]; // свободные члены уравнений системы

        // Устанавливаем для каждого пикселя свободный член b (на основе градиента вставляемого изображения и граничных пикселей исходного)
        for (int i = 0; i < pixel.size(); i++) {

            lines[i] = new MatrixLine();

            Pixel pix = pixel.get(i);
            int x = pix.x;
            int y = pix.y;
            int b_ = 0;

            int N = 0;

            b_ += grad[x + y * w];
            for (int p = 0; p < D.length; p++) {
                int x_ = x + D[p][0];
                int y_ = y + D[p][1];
                if (x_ < 0 || y_ < 0 || x_ >= w || y_ >= h) {
                    continue;
                }
                int ind = x_ + y_ * w;
                if (alpha[ind] > 0) {
                    if (h_alpha[ind] > 0) {
                        b_ += h_color[ind];
                    } else {
                        lines[i].add(index_to_number[ind], -1);
                    }
                    N++;
                }
            }

            b[i] = b_;

            lines[i].add(i, N);
        }

        double[] r = solve(lines, b, lines.length); // решаем систему методом Гаусса-Зейделя
        for (int i = 0; i < r.length; i++) {
            r[i] /= 1.2;
            res[number_to_index[i]] = (int) Math.abs(r[i]) % 256; // преобразуем цвета
        }

        return res;
    }

    private class Pixel {

        int x, y;

        public Pixel(int a, int b) {
            x = a;
            y = b;
        }
    }

    private int[] getGradient(int[] pix, int w) {
        int N = 3;
        int h = pix.length / w;
        int[] pq = new int[pix.length];
        int[][] CONV = {{0, -1, 0}, {-1, 4, - 1}, {0, -1, 0}};
        for (int j = 0; j < h; j++) {
            int pqa = j * w;
            for (int i = 0; i < w; i++) {
                int px = 0;//, py = 0;
                for (int x = -N / 2; x < N / 2 + 1; x++) {
                    for (int y = -N / 2; y < N / 2 + 1; y++) {
                        int x_r = Math.abs(x + i);
                        int y_r = Math.abs(y + j);
                        if (x_r >= w) {
                            x_r = w - 1;
                        }
                        if (y_r >= h) {
                            y_r = h - 1;
                        }
                        px += pix[x_r + y_r * w] * CONV[x + N / 2][y + N / 2];
                        //  py += pix[x_r + y_r * w] * SOBEL[y + N / 2][x + N / 2];
                    }
                }
               // px /= 8;

                // py /= 8;
                // px = (int) Math.sqrt(px * px + py * py);
                pq[i + pqa] = px;
            }
        }
        return pq;
    }

    private Image getPoisson(Image img, Image fon) {
        int[] pix = getPixels(img);

        // pix = getGrayScale(pix);
        // массив пикселей изображения размера img.getWidth() x img.getHeight() 
        // со значениями пикселов фонового изображения img на границах
        //pix = mergeChannels(getChannel(pix, 24), getGradient(getChannel(pix, 16), img.getWidth(null)), getGradient(getChannel(pix, 8), img.getWidth(null)), getGradient(getChannel(pix, 0), img.getWidth(null)));
        int[] h = getH(fon, getChannel(pix, 24), img.getWidth(null), ANCHOR_X, ANCHOR_Y);
        int[] pixa = getChannel(pix, 24);
        int[] ha = getChannel(h, 24);
        pix = mergeChannels(pixa, runPoisson(pixa, getChannel(pix, 16), ha, getChannel(h, 16), img.getWidth(null)), runPoisson(pixa, getChannel(pix, 8), ha, getChannel(h, 8), img.getWidth(null)), runPoisson(pixa, getChannel(pix, 0), ha, getChannel(h, 0), img.getWidth(null)));
        return createImage(new MemoryImageSource(img.getWidth(null), img.getHeight(null), pix, 0, img.getWidth(null)));
    }

    private double[] solve(MatrixLine[] a, double[] b, int n) {
        double[] x = new double[n];
        double[] p = new double[n];
        int i1 = 0;
        do {
            System.arraycopy(x, 0, p, 0, n);

            for (int i = 0; i < n; i++) {
                double var = 0;
                int[][] d = a[i].data;
                int ind_m = 0;
                for (int m = 0; m < a[i].count; m++) {
                    int j = d[m][0];
                    if (j < i) {
                        var += (d[m][1] * x[j]);
                    } else if (j > i) {
                        var += (d[m][1] * p[j]);
                    } else {
                        ind_m = m;
                    }
                }
                x[i] = (b[i] - var) / d[ind_m][1];
            }
            i1++;
        } while (!converge(x, p));
        System.out.println("Итераций " + i1);
        return x;
    }

    boolean converge(double[] xk, double[] xkp) {
        double norm = 0;
        for (int i = 0; i < xk.length; i++) {
            norm += (xk[i] - xkp[i]) * (xk[i] - xkp[i]);
        }
        return Math.sqrt(norm) < eps;
    }

    private class MatrixLine {

        public int[][] data;
        public int count;

        public MatrixLine() {
            data = new int[5][2];
        }

        public void add(int x, int el) {
            data[count][0] = x;
            data[count++][1] = el;
        }
    }

}
