# Print an rmst_power result

Prints the model metadata and power table returned by
[`rmst.power()`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/rmst.power.md).

Returns a summary object for printing and further inspection.

Prints and returns the stored `ggplot2` object.

## Usage

``` r
# S3 method for class 'rmst_power'
print(x, ...)

# S3 method for class 'rmst_power'
summary(object, ...)

# S3 method for class 'rmst_power'
plot(x, ...)
```

## Arguments

- x:

  An object returned by
  [`rmst.power()`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/rmst.power.md).

- ...:

  Additional arguments passed to or from other methods.

- object:

  An object returned by
  [`rmst.power()`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/rmst.power.md).

## Value

The input object `x`, invisibly.

An object of class `"summary.rmst_power"`.

A `ggplot2` object, invisibly.
