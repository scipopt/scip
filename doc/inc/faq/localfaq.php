<!-- this file is used to generate the local doxygen documentation -->
<!-- using make doc from within scip or soplex -->

<style>
@media (prefers-color-scheme: light) {
.reveal:hover {
    text-shadow: 1px 1px 1px #777;
}
.answer {
    background-color: #fff;
}
}
.answer {
    padding-left:   1em;
}
</style>

<?php include('faqdata.php'); ?>
<?php
//Output a table of contents

foreach ($faq as $section) {
  echo '<h3><a class="reveal_faq" href="#'.$section['label'].'"><span class="fa fa-caret-right"></span> '.$section['title'].'</a></h3>';
  echo '<ol>';

    foreach($section['content'] as $item) {
      $label = $item['label'];
  ?>
  <li>
    <a class="reveal_faq" href="#<?php echo $label ?>">
      <?php echo $item['question'] ?>
    </a>
  </li>
  <?php
  }

  echo "</ol>\n";
}
echo "<hr />";
?>

<?php
// Content: questions and answers

foreach ($faq as $section) {
  echo '<h3 id="'.$section['label'].'" class="anchor"><span class="fa fa-caret-right"></span> '.$section['title'].'<a href="#" class="pull-right"><span title="go to top" class="fa fa-caret-up"></span></a></h3>';
  echo '<ol>';

    foreach($section['content'] as $item) {
      $label = $item['label'];
  ?>
  <li id="<?php echo $label ?>" class="anchor">
    <h4>
      <a class="reveal_faq" href="#<?php echo $label ?>"><?php echo $item['question'] ?></a>

      <a href="#" class="pull-right"><span class="fa fa-caret-up" title="go to top"></span></a>
    </h4>
    <div id="<?php echo $label ?>_ans" class="answer">
      <?php echo $item['answer'] ?>
    </div>
  </li>
  <?php
  }

  echo "</ol>\n";
}
?>
<!--WE don't use this stuff because it is uncontrollable
<script>
$(".reveal").click(function() {
    var questionId = "#"+this.id+"_ans";
    $(questionId).toggle("fast");
});
</script>-->
