# Run AFTER lysozyme_analysis.py (needs results/residue_table.csv)
# install.packages("ggplot2")   # only once

library(ggplot2)

df <- read.csv("results/residue_table.csv")

p <- ggplot(df, aes(x = contacts, y = rsasa, colour = location)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE,
              colour = "black", linewidth = 0.6) +
  geom_point(data = subset(df, catalytic_site == "True" | catalytic_site == TRUE),
             shape = 21, size = 3.5, colour = "red", stroke = 1) +
  scale_colour_manual(values = c(buried = "#1f77b4", exposed = "#ff7f0e")) +
  labs(title = "Lysozyme (1LYZ): burial vs solvent accessibility",
       subtitle = "Red rings = residues near the catalytic site",
       x = "Contacts within 8 A (CA-CA)",
       y = "Relative SASA",
       colour = "Location") +
  theme_minimal(base_size = 12)

ggsave("results/rsasa_vs_contacts_R.png", p, width = 6.5, height = 5, dpi = 300)